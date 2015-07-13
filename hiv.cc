/*
 * Main PhD code for comparing the effect of model changes on outputs.
 *

 ### Model A

 Deterministic population growth E1
 Deterministic population deaths by HIV and ARV status E2
     (5-stage Granich model)
 Deterministic HIV infection (including by risk strata) 5
 Including:
 - Deterministic disease progression
 - Deterministic go onto ARVs

 ### Model B

 Deterministic population immigration
 Deterministic population emigration
 Deterministic population births
 Deterministic population deaths by HIV and ARV status (5-stage Granich model)
 Stochastic HIV infection 6

 ### Model C

 Deterministic population immigration
 Deterministic population emigration
 Deterministic population births
 Stochastic population deaths by HIV and ARV status (infection time) 7
 Deterministic HIV infection

 ### Model D

 Deterministic population immigration
 Deterministic population emigration
 Deterministic population births
 Stochastic population deaths by HIV and ARV status (infection time)
 Stochastic HIV infection

 ### Model E

 Stochastic population immigration 8
 Stochastic population emigration 9
 Stochastic population births 10
 Stochastic population deaths by HIV and ARV status (infection time)
 Stochastic HIV infection

 ### Model F

 Stochastic population immigration
 Stochastic population emigration
 Stochastic population births
 Stochastic population deaths by HIV and ARV status (infection time)
 Stochastic HIV infection using Leigh's method 11

 ### Model G

 Stochastic population immigration
 Stochastic population emigration
 Stochastic population births
 Stochastic population deaths by HIV and ARV status (infection time)
 Stochastic random matching for HIV infection 12


 ### Model H

 Stochastic population immigration
 Stochastic population emigration
 Stochastic population births
 Stochastic population deaths by HIV and ARV status (infection time)
 Stochastic HIV infection using continuous heterogeneity k-random 13

 ### Model I

 Stochastic population immigration
 Stochastic population emigration
 Stochastic population births
 Stochastic population deaths by HIV and ARV status (infection time)
 Stochastic HIV infection using weighted shuffle 14


 ### Model J

 Stochastic population immigration
 Stochastic population emigration
 Stochastic population births
 Stochastic population deaths by HIV and ARV status (infection time)
 Stochastic HIV infection using cluster shuffle 15

*/

#include "sim.hh"

void calcVariablesEvent(sim::Simulation &s)
{
  size_t num_susceptible = 0, num_infected = 0;
  size_t num_stage_1 = 0, num_stage_2 = 0, num_stage_3 = 0, num_stage_4 = 0;
  size_t num_arvs = 0;

  for (auto & a: s.agents) {
    HIVAgent *agent = (HIVAgent *) a;
    if(agent->hiv == 0) {
      ++num_susceptible;
    } else {
      ++num_infected;
      switch(agent->hiv) {
      case 1:
	++num_stage_1;
	break;
      case 2:
	++num_stage_2;
	break;
      case 3:
	++num_stage_3;
	break;
      case 4:
	++num_stage_4;
	break;
      case 5:
	++num_arvs;
	break;
      default:
	throw InvalidData(std::string("Invalid HIV value."));
      }
    }
  }
  s.context.set("_AGENTS",s.agents.size());
  s.context.set("_AGENTS_0", num_susceptible);
  s.context.set("_AGENTS_INFECTED", num_infected);
  s.context.set("_AGENTS_1", num_stage_1);
  s.context.set("_AGENTS_2", num_stage_2);
  s.context.set("_AGENTS_3", num_stage_3);
  s.context.set("_AGENTS_4", num_stage_4);
  s.context.set("_AGENTS_5", num_arvs);
}

/* E1 */

void increasePopulationDeterministic(sim::Simulation &s)
{
  std::uniform_real_distribution<double> uni;
  unsigned num_agents_to_create;
  double growth_rate = s.context("GROWTH");
  double initial_age = s.context("INITIAL_AGE");
  s.context.set_if_not_set("_PARTIAL_GROWTH", { 0.0 });
  double &partial_growth = s.context.get("_PARTIAL_GROWTH")[0];
  double growth;
  size_t num_agents = s.context("_AGENTS");

  growth = std::max(0.0, growth_rate * (num_agents + partial_growth));
  partial_growth = growth - std::floor(growth);

  if (growth > num_agents) {
    num_agents_to_create = std::floor(growth) - num_agents;

    for (unsigned i = 0; i < num_agents_to_create; ++i) {
      HIVAgent *a = create_hiv_agent(s.context);
      a->sex = (i % 2 == 0) ? sim::MALE : sim::FEMALE;
      a->dob = s.current_date - initial_age;
      a->hiv = 0;
      s.agents.push_back(a);
    }
  }
}

/*
 * E2
 *
 * Stage 0: HIV-negative
 * Stage 1 to 4: WHO stages
 * Stage 5: ARVs
 */

void decreasePopulationDeterministic(sim::Simulation &s)
{
  for (int stage = 0; stage < 6; ++stage) {
    unsigned num_agents_to_remove = 0;
    std::string parm = "DECLINE_RATE_" + std::to_string(stage);
    double decline_rate = 1.0 - s.context(parm.c_str());
    parm = "_PARTIAL_DECLINE_" + std::to_string(stage);
    s.context.set_if_not_set(parm.c_str(), { 0.0 });
    double &partial_decline = s.context.get(parm.c_str())[0];
    double decline;
    parm = "_AGENTS_" + std::to_string(stage);
    size_t num_agents = s.context(parm.c_str());

    decline = std::max(0.0, decline_rate * (num_agents + partial_decline));
    partial_decline = decline - std::floor(decline);

    if (decline < num_agents)
      num_agents_to_remove = num_agents - std::floor(decline);

    unsigned i = 0;
    for (auto it = s.agents.begin();
	 it != s.agents.end() && i < num_agents_to_remove;
	 ++it, ++i) {
      HIVAgent *agent = (HIVAgent *) *it;
      if ( agent->hiv == stage)
	agent->die(s, "DETERMINISTIC");
    }
  }
  s.remove_dead_agents();
}


void emigrationDiffEqEvent(sim::Simulation &s)
{
  std::uniform_real_distribution<double> uni;
  unsigned num_agents_to_remove;
  double decline_rate = s.context("DECLINE");
  double stdev = s.context("DECLINE_STDEV");
  std::normal_distribution<double> normal(decline_rate, stdev);
  s.context.set_if_not_set("_PARTIAL_DECLINE", { 0.0 });
  double &partial_decline = s.context.get("_PARTIAL_DECLINE")[0];
  double decline, floor_decline;
  size_t num_agents = s.context("_AGENTS");

  decline_rate = normal(sim::rng);
  decline = decline_rate * (num_agents + partial_decline);

  if (fabs(round(decline) - decline) < 0.000000000001) {
    if (decline >= round(decline))
      floor_decline = floor(decline);
    else
      floor_decline = ceil(decline);
  } else {
    floor_decline = floor(decline);
  }

  partial_decline = decline - floor_decline;

  if (num_agents < floor_decline)
    num_agents_to_remove = num_agents;
  else
    num_agents_to_remove = num_agents -  (size_t) floor_decline;

  for (unsigned i = 0; i < num_agents_to_remove; ++i) {
    s.aged_out_agents.push_back(s.agents.back());
    s.agents.pop_back();
  }
}

void report(sim::Simulation &s)
{
  unsigned num_alive = s.agents.size();
  unsigned num_susceptible = count_if(s.agents.begin(), s.agents.end(),
				      [](const sim::Agent *a) {
					const HIVAgent *b = (HIVAgent *) a;
					return b->hiv == 0;
				      });
  unsigned num_infected = count_if(s.agents.begin(), s.agents.end(),
				   [](const sim::Agent *a) {
				     const HIVAgent *b = (HIVAgent *) a;
				     return b->hiv > 0;
				   });
  unsigned num_arvs = count_if(s.agents.begin(), s.agents.end(),
			       [](const sim::Agent *a) {
				 const HIVAgent *b = (HIVAgent *) a;
				 return b->hiv == 5;
			       });
  unsigned num_dead = s.dead_agents.size();
  unsigned num_moved_out = s.aged_out_agents.size();

  std::stringstream ss;
  ss << s.simulation_num << ", " << s.current_date << ", alive, "
     << num_alive << std::endl;
  ss << s.simulation_num << ", " << s.current_date << ", susceptible, "
     << num_susceptible << std::endl;
  ss << s.simulation_num << ", " << s.current_date << ", infected, "
     << num_infected << std::endl;
  ss << s.simulation_num << ", " << s.current_date << ", on_arvs, "
     << num_arvs << std::endl;
  ss << s.simulation_num << ", " << s.current_date << ", dead, "
     << num_dead << std::endl;
  ss << s.simulation_num << ", " << s.current_date << ", moved_out, "
     << num_moved_out << std::endl;

  std::cout << ss.str();
}

void testSti(tst::TestSeries &t)
{

  sim::Simulation(sim::Options()
		  .events({sim::advanceTimeEvent
			, calcVariablesEvent
			, increasePopulationDeterministic})
		  .afterEachSimulation([&t](sim::Simulation &s) {
		      TESTEQ(t, s.agents.size(), 6727, "Population growth");
		    })
		  .agentCreate(create_hiv_agent)
		  .timeAdjust("GROWTH", 1.0, 0, 1, sim::COMPOUND)
		  .parameter("NUM_AGENTS", {1000.0} )
		  .parameter("TIME_STEP", {1.0 / 365.25} )
		  .parameter("INITIAL_AGE", {15.0} )
		  .parameter("GROWTH", {1.1} )
		  .parameter("GROWTH_STDEV", {0.0} )
		  .parameter("RISK_SUSCEPTIBLE", {0.0} )
		  .parameter("RISK_SUSCEPTIBLE_STDEV", {0.0} )
		  .parameter("RISK_UNINFECTIOUS", {0.0} )
		  .parameter("RISK_UNINFECTIOUS_STDEV", {0.0} )
		  .parameter("PARTNER_FORMATION_RATE", {0.0} )
		  .parameter("PARTNER_FORMATION_RATE_STDEV", {0.0} )
		  .parameter("PREVALENCE", {0.0} )).simulate();

  sim::Simulation(sim::Options()
		  .events({sim::advanceTimeEvent
			, calcVariablesEvent
			, decreasePopulationDeterministic})
		  .afterEachSimulation([&t](sim::Simulation &s) {
		      TESTEQ(t, s.agents.size(), 817, "Population decline");
		    })
		  .agentCreate(create_hiv_agent)
		  .timeAdjust("DECLINE_RATE_0", 1.0, 0, 1, sim::PROBABILITY)
		  .timeAdjust("DECLINE_RATE_1", 1.0, 0, 1, sim::COMPOUND)
		  .timeAdjust("DECLINE_RATE_2", 1.0, 0, 1, sim::COMPOUND)
		  .timeAdjust("DECLINE_RATE_3", 1.0, 0, 1, sim::COMPOUND)
		  .timeAdjust("DECLINE_RATE_4", 1.0, 0, 1, sim::COMPOUND)
		  .timeAdjust("DECLINE_RATE_5", 1.0, 0, 1, sim::COMPOUND)
		  .parameter("NUM_AGENTS", {1000.0} )
		  .parameter("DECLINE_RATE_0", {0.01} )
		  .parameter("DECLINE_RATE_1", {0.02} )
		  .parameter("DECLINE_RATE_2", {0.03} )
		  .parameter("DECLINE_RATE_3", {0.1} )
		  .parameter("DECLINE_RATE_4", {0.5} )
		  .parameter("DECLINE_RATE_5", {0.015} )
		  .parameter("TIME_STEP", {1.0 / 365.25 } )
		  .parameter("INITIAL_AGE", {15.0} )
		  .parameter("GROWTH", {1.1} )
		  .parameter("HIV_PREVALENCE", {0.0} )).simulate();


  sim::Simulation(sim::Options()
		  .events({sim::advanceTimeEvent
			, calcVariablesEvent
			, emigrationDiffEqEvent})
		  .afterEachSimulation([&t](sim::Simulation &s) {
		      TESTEQ(t, s.agents.size(), 423, "Population decline");
		    })
		  .agentCreate(create_hiv_agent)
		  .timeAdjust("DECLINE", 1.0, 0, 1, sim::COMPOUND)
		  .parameter("NUM_AGENTS", {10000.0} )
		  .parameter("TIME_STEP", {1.0 / 365.25} )
		  .parameter("NUM_YEARS", {30.0} )
		  .parameter("INITIAL_AGE", {15.0} )
		  .parameter("DECLINE", {0.9} )
		  .parameter("DECLINE_STDEV", {0.0} )
		  .parameter("RISK_SUSCEPTIBLE", {0.0} )
		  .parameter("RISK_SUSCEPTIBLE_STDEV", {0.0} )
		  .parameter("RISK_UNINFECTIOUS", {0.0} )
		  .parameter("RISK_UNINFECTIOUS_STDEV", {0.0} )
		  .parameter("PARTNER_FORMATION_RATE", {0.0} )
		  .parameter("PARTNER_FORMATION_RATE_STDEV", {0.0} )
		  .parameter("HIV_PREVALENCE", {0.0} )).simulate();



}


int main(int argc, char **argv)
{

  // Model A

  sim::Simulation(sim::Options()
		  .additionalTests({testSti})
		  .beforeEachSimulation(report)
		  .events({sim::advanceTimeEvent
			, calcVariablesEvent
			// , increasePopulationDeterministic
			// , report
			})
		  .afterEachSimulation(report)
		  .commandLine(argc, argv)
		  .agentCreate(create_hiv_agent)).simulate();
  return 0;
}
