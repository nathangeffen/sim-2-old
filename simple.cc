#include "sim.hh"

sim::HIVAgent *create_hiv_agent_simple(sim::Context &c)
{
  std::uniform_real_distribution<double> uni;
  double hiv_infection_time;
  sim::HIVAgent *a = new sim::HIVAgent(c);
  a->sex = uni(sim::rng) < c("PROB_MALE") ?
			   sim::MALE : sim::FEMALE;

  std::uniform_real_distribution<double>
		      uni_age(c("EARLIEST_BIRTH_DATE"),
			      c("LATEST_BIRTH_DATE"));
  a->dob = uni_age(sim::rng);

  // HIV negative CD4 count
  std::normal_distribution<double> cd4_distribution(c("CD4_HIV_NEG_AVG"),
						    c("CD4_HIV_NEG_STDEV"));
  a->cd4 = cd4_distribution(sim::rng);
  a->cd4 = std::max(c("CD4_MIN"), a->cd4);
  a->cd4 = std::min(c("CD4_MAX"), a->cd4);
  a->hiv = 0;
  if (uni(sim::rng) < c("HIV_PREVALENCE")) { // HIV positive
    a->hiv = 1.0;
    std::uniform_real_distribution<double>
      uni_hiv_time(0, c("MAX_INFECTION_PERIOD"));
    hiv_infection_time = uni_hiv_time(sim::rng);
    a->hiv_infection_date = c("START_DATE") - hiv_infection_time;
    if (a->hiv_infection_date < a->dob) {
      a->hiv_infection_date = a->dob;
      hiv_infection_time = c("START_DATE") - a->dob;
    }
    a->cd4 -= std::max(0.0,
		       hiv_infection_time * c("AVG_CD4_DECLINE_UNCORRECTED"));
  }
  a->riskiness = uni(sim::rng);
  a->orientation = 1.0;
  a->num_partners = 0;

  if (uni(sim::rng) < c("PROB_CIRCUMCISED"))
    a->circumcised = true;

  return a;
}

void cd4Event(sim::Simulation &s)
{
  static thread_local double avg_cd4_decline = s.context("AVG_CD4_DECLINE");
  for (auto a : s.agents) {
    sim::HIVAgent *agent = (sim::HIVAgent *) a;
    if (agent->hiv)
      agent->cd4 = std::max(0.0, agent->cd4 - avg_cd4_decline);
  }
}

/* Simple SIR model for HIV infection.*/
void hivEvent(sim::Simulation &s)
{
  double incidence = s.context("HIV_INCIDENCE");
  double risk_infection, S, I;
  unsigned number_hiv = 0, number_no_hiv;

  for (auto a : s.agents) {
    sim::HIVAgent *agent = (sim::HIVAgent *)  a;
    if (agent->hiv)
      ++number_hiv;
  }

  number_no_hiv = s.agents.size() - number_hiv;
  S = (double) number_no_hiv / s.agents.size();
  I = (double) number_hiv / s.agents.size();
  risk_infection = incidence * S * I;

  for (auto a : s.agents) {
    sim::HIVAgent *agent = (sim::HIVAgent *) a;
    if (agent->hiv == 0) {
      if (sim::is_event(risk_infection)) {
	agent->hiv = 1;
	agent->hiv_infection_date = s.current_date;
      }
    }
  }
}

void deathEvent(sim::Simulation &s)
{
  double risk_of_death, hiv_risk_of_death;
  unsigned age;
  std::vector<double> & male_mort =
    s.context.get("RISK_DEATH_MALE");
  std::vector<double> & female_mort=
    s.context.get("RISK_DEATH_FEMALE");
  std::vector<double> & hiv_mort =
    s.context.get("HIV_RISK_DEATH");

  for (auto a : s.agents) {
    sim::HIVAgent *agent = (sim::HIVAgent *) a;
    age = agent->age(s);
    // Risk of death for HIV- due to age
    if (agent->sex == sim::MALE)
      risk_of_death = male_mort[age];
    else
      risk_of_death = female_mort[age];
    if (agent->hiv == 0) {
      if (sim::is_event(risk_of_death))
	agent->die(s, "AGE");
    } else  {
    // HIV risk of death
      size_t index;
      if (age < 5.0)
	index = 0;
      else if (age < 15.0)
	index = 6;
      else if (age < 25.0)
	index = 12;
      else if (age < 35.0)
	index = 18;
      else if (age < 45.0)
	index = 24;
      else if (age < 55.0)
	index = 30;
      else
	index = 36;
      if (agent->cd4 < 50.0)
	;
      else if (agent->cd4 < 100.0)
	++index;
      else if (agent->cd4 < 200.0)
	index += 2;
      else if (agent->cd4 < 350.0)
	index += 3;
      else if (agent->cd4 < 500.0)
	index += 4;
      else
	index += 5;
      hiv_risk_of_death = std::max(hiv_mort[index], risk_of_death);
      if (sim::is_event(hiv_risk_of_death)) {
	agent->die(s, "HIV");
      }
    }
  }
  s.remove_dead_agents();
}

void print_agents(sim::Simulation &s)
{
  std::cout << "ALIVE" << std::endl;
  for (auto a : s.agents) {
    sim::HIVAgent *agent = (sim::HIVAgent *) a;
    std::cout << *agent << std::endl;
  }
  std::cout << "DEAD" << std::endl;
  for (auto a : s.dead_agents) {
    sim::HIVAgent *agent = (sim::HIVAgent *) a;
    std::cout << *agent << std::endl;
  }
}

double calcPrevalence(sim::Simulation &s)
{
  unsigned num_hiv_pos = 0;

  for (auto a : s.agents) {
    sim::HIVAgent *agent = (sim::HIVAgent *) a;
    if (agent->hiv)
      ++num_hiv_pos;
  }
  return (double) num_hiv_pos / s.agents.size();
}

void interimReport(sim::Simulation &s)
{
  std::stringstream ss;
  ss << s.simulation_num << ", " << s.current_iteration
     << ", HIV prevalence, " << calcPrevalence(s) << std::endl;
  std::cout << ss.str();
}

void report(sim::Simulation &s)
{
  unsigned num_alive_hiv = 0;
  unsigned num_died_hiv = 0;
  unsigned num_dead_hiv = 0;
  double avg_hiv_years = 0.0;
  double avg_hiv_years_alive = 0.0;
  double avg_age_all_alive = 0.0;
  double avg_age_hiv_negative_alive = 0.0;
  double avg_age_hiv_positive_alive = 0.0;
  double avg_age_all_dead = 0.0;
  double avg_age_hiv_negative_dead = 0.0;
  double avg_age_hiv_positive_dead = 0.0;
  double avg_cd4_hiv_positive_alive = 0.0;
  double avg_cd4_hiv_positive_dead = 0.0;

  double prop_dead = (double) s.dead_agents.size() /
    (s.agents.size() + s.dead_agents.size());
  std::stringstream ss;
  ss << s.simulation_num << ", alive, " << s.agents.size() << std::endl;
  ss << s.simulation_num << ", dead, "
     << s.dead_agents.size() << std::endl;
  ss << s.simulation_num << ", proportion dead, "
     << prop_dead << std::endl;

  for (auto a : s.agents) {
    sim::HIVAgent *agent = (sim::HIVAgent *) a;
    avg_age_all_alive += s.current_date - agent->dob;
    if (agent->hiv) {
      ++num_alive_hiv;
      avg_age_hiv_positive_alive += s.current_date - agent->dob;
      avg_hiv_years_alive += s.current_date - agent->hiv_infection_date;
      avg_cd4_hiv_positive_alive += agent->cd4;
    } else {
      avg_age_hiv_negative_alive += s.current_date - agent->dob;
    }
  }
  avg_hiv_years_alive /= num_alive_hiv;
  avg_age_hiv_positive_alive /= num_alive_hiv;
  avg_age_hiv_negative_alive /= (s.agents.size() - num_alive_hiv);
  avg_age_all_alive /= s.agents.size();
  avg_cd4_hiv_positive_alive /= num_alive_hiv;

  for (auto a : s.dead_agents) {
    sim::HIVAgent *agent = (sim::HIVAgent *) a;
    avg_age_all_dead += agent->dod - agent->dob;
    if (agent->hiv) {
      ++num_dead_hiv;
      avg_age_hiv_positive_dead += agent->dod - agent->dob;
      avg_hiv_years += agent->dod - agent->hiv_infection_date;
      avg_cd4_hiv_positive_dead += agent->cd4;
      if (agent->cause == "HIV")
	++num_died_hiv;
    } else {
      avg_age_hiv_negative_dead += agent->dod - agent->dob;
    }
  }
  avg_hiv_years /= num_dead_hiv;
  avg_age_hiv_positive_dead /= num_dead_hiv;
  avg_age_hiv_negative_dead /= (s.dead_agents.size() - num_dead_hiv);
  avg_age_all_dead /= s.dead_agents.size();
  avg_cd4_hiv_positive_dead /= num_dead_hiv;

  ss << s.simulation_num << ", HIV alive, " << num_alive_hiv << std::endl;
  ss << s.simulation_num << ", avg years alive with HIV, "
     << avg_hiv_years_alive << std::endl;
  ss << s.simulation_num << ", HIV dead, " << num_dead_hiv << std::endl;
  ss << s.simulation_num << ", avg years till death if HIV, "
     << avg_hiv_years << std::endl;
  ss << s.simulation_num << ", HIV dead who died of HIV, "
     << num_died_hiv << std::endl;
  ss << s.simulation_num << ", avg age alive, "
     << avg_age_all_alive << std::endl;
  ss << s.simulation_num << ", avg age dead, "
     << avg_age_all_dead << std::endl;
  ss << s.simulation_num << ", avg HIV positive age alive, "
     << avg_age_hiv_positive_alive << std::endl;
  ss << s.simulation_num << ", avg HIV negative age alive, "
     << avg_age_hiv_negative_alive << std::endl;
  ss << s.simulation_num << ", avg HIV positive age dead, "
     << avg_age_hiv_positive_dead << std::endl;
  ss << s.simulation_num << ", avg HIV negative age dead, "
     << avg_age_hiv_negative_dead << std::endl;
  ss << s.simulation_num << ", avg HIV CD4 alive, "
     << avg_cd4_hiv_positive_alive << std::endl;
  ss << s.simulation_num << ", avg HIV CD4 dead, "
     << avg_cd4_hiv_positive_dead << std::endl;

  std::cout << ss.str();
}


int main(int argc, char **argv)
{
  sim::Simulation(sim::Options()
		  .events({sim::advanceTimeEvent
			, cd4Event
			, hivEvent
			, deathEvent
			, interimReport})
		  .afterEachSimulation(report)
		  .commandLine(argc, argv)
		  .agentCreate(create_hiv_agent_simple)).simulate();
  return 0;
}
