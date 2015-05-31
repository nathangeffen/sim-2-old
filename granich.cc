#include "sim.hh"

sim::HIVAgent *createHivAgentGranichInit(sim::Context &c)
{
  std::uniform_real_distribution<double> uni;
  sim::HIVAgent *a = new sim::HIVAgent(c);

  a->sex = uni(sim::rng) < c("PROB_MALE") ?
			   sim::MALE : sim::FEMALE;
  { // dob
    std::uniform_real_distribution<double>
      uni_age(c("EARLIEST_BIRTH_DATE"),
	      c("LATEST_BIRTH_DATE"));
    a->dob = uni_age(sim::rng);
  }
  a->cd4 = 1000;
  a->hiv = 0;
  return a;
}

void createHivAgentsGranich(sim::Simulation &s)
{
  std::uniform_real_distribution<double> uni;
  unsigned num_agents_to_create;
  double beta = s.context("BETA");
  double stdev = s.context("BETA_STDEV");
  std::normal_distribution<double> normal(beta, stdev);
  double new_beta = normal(sim::rng);
  double initial_age = s.context("INITIAL_AGE");

  num_agents_to_create = round(new_beta * s.agents.size());
  std::cout << "D0: " << new_beta << " " << num_agents_to_create << std::endl;

  for (unsigned i = 0; i < num_agents_to_create; ++i) {
    sim::HIVAgent *a = new sim::HIVAgent(s.context);
    a->sex = uni(sim::rng) < s.context("PROB_MALE") ?
			     sim::MALE : sim::FEMALE;
    a->dob = s.current_date - initial_age;
    a->cd4 = 1000.0;
    a->hiv = 0;
    s.agents.push_back(a);
  }
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
			, createHivAgentsGranich})
		  .afterEachSimulation(report)
		  .commandLine(argc, argv)
		  .agentCreate(createHivAgentGranichInit)).simulate();
  return 0;
}
