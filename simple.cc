#include "sim.hh"

void cd4Event(sim::Simulation &s)
{
  static thread_local double avg_cd4_decline = s.common("AVG_CD4_DECLINE");
  for (auto agent : s.agents) {
    if (agent->hiv) {
      agent->cd4 = std::max(0.0,
			    agent->cd4 -
			    avg_cd4_decline * s.time_step);
    }
  }
}


void deathEvent(sim::Simulation &s)
{
  static thread_local std::vector<double> & male_mort =
    s.common.get("RISK_DEATH_MALE");
  static thread_local std::vector<double> & female_mort=
    s.common.get("RISK_DEATH_FEMALE");
  static thread_local std::vector<double> & hiv_mort =
    s.common.get("HIV_RISK_DEATH");
  static thread_local size_t size = hiv_mort.size();

  for (auto agent : s.agents) {
    unsigned age = agent->age(s);
    double risk_of_death;

    // Risk of death for everyone due to age
    if (agent->sex == sim::MALE)
      risk_of_death = male_mort[age];
    else
      risk_of_death = female_mort[age];
    if (sim::is_event(risk_of_death))
      agent->die(s, "AGE");

    // HIV risk of death
    if (agent->alive && agent->hiv) {
      for (size_t i = size; i > 0; i -= 2) {
	if (agent->cd4 < hiv_mort[i - 2]) {
	  double hiv_risk_of_death = hiv_mort[i - 1];
	  if (sim::is_event(hiv_risk_of_death)) {
	    agent->die(s, "HIV");
	    break;
	  }
	}
      }
    }
  }
  s.remove_dead_agents();
}


void report(sim::Simulation &s)
{
  double prop_dead = (double) s.dead_agents.size() /
    (s.agents.size() + s.dead_agents.size());
  std::cout << s.simulation_num << ", alive, " << s.agents.size() << std::endl;
  std::cout << s.simulation_num << ", dead, "
	    << s.dead_agents.size() << std::endl;
  std::cout << s.simulation_num << ", proportion dead, "
	    << prop_dead << std::endl;
}

int main(int argc, char *argv[])
{
  sim::Common c;
  sim::Simulation s(c);

  s.simulate({sim::advanceTimeEvent, cd4Event, deathEvent}, report,
	     argc, argv);

  return 0;
}
