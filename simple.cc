#include "sim.hh"

void cd4Event(sim::Simulation &s)
{
  static thread_local double avg_cd4_decline = s.context("AVG_CD4_DECLINE");
  for (auto agent : s.agents)
    if (agent->hiv)
      agent->cd4 = std::max(0.0, agent->cd4 - avg_cd4_decline);
}

void hivEvent(sim::Simulation &s)
{

  return;
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

  for (auto agent : s.agents) {
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
	break;
      }
    }
  }
  s.remove_dead_agents();
}

void report(sim::Simulation &s)
{
  double prop_dead = (double) s.dead_agents.size() /
    (s.agents.size() + s.dead_agents.size());
  std::stringstream ss;
  ss << s.simulation_num << ", alive, " << s.agents.size() << std::endl;
  ss << s.simulation_num << ", dead, "
     << s.dead_agents.size() << std::endl;
  ss << s.simulation_num << ", proportion dead, "
     << prop_dead << std::endl;
  std::cout << ss.str();
}

int main(int argc, char *argv[])
{
  sim::Simulation simulation;

  simulation.simulate({
      sim::advanceTimeEvent,
	cd4Event,
	hivEvent,
	deathEvent},
    report, argc, argv);

  return 0;
}
