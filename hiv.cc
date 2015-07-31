/*
 * Main PhD code for comparing the effect of model changes on outputs.
 *

 ### Model A

 Deterministic population growth E1
 Deterministic population deaths by HIV and ARV status E2
     (5-stage Granich model)
 Deterministic disease progression E3
 Deterministic HIV infection (including by risk strata) 5 E4
 Deterministic ARVs E5

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
  size_t num_susceptible = 0, num_susceptible_high = 0, num_susceptible_low = 0;
  size_t num_infected = 0;
  size_t num_stage_1 = 0, num_stage_2 = 0, num_stage_3 = 0, num_stage_4 = 0;
  size_t num_arvs = 0;

  for (auto & a: s.agents) {
    HIVAgent *agent = (HIVAgent *) a;
    if(agent->hiv == 0) {
      ++num_susceptible;
      if (agent->risk == 0)
	++num_susceptible_low;
      else
	++num_susceptible_high;
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
  s.context.set("_AGENTS_0_LOW", num_susceptible_low);
  s.context.set("_AGENTS_0_HIGH", num_susceptible_high);
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
  unsigned num_agents_to_create = 0;
  double growth_rate = 1.0 + s.context("GROWTH");
  double initial_age = s.context("INITIAL_AGE");
  s.context.set_if_not_set("_PARTIAL_GROWTH", { 0.0 });
  double &partial_growth = s.context.get("_PARTIAL_GROWTH")[0];
  double growth;
  size_t num_agents = s.context("_AGENTS");

  growth = std::max(0.0, growth_rate * (num_agents + partial_growth));
  partial_growth = growth - std::floor(growth);

  if (growth > num_agents) {
    num_agents_to_create = std::floor(growth) - num_agents;
    double high_risk = s.context("HIGH_RISK_PROPORTION");
    size_t num_high_risk = high_risk * (double) num_agents_to_create;
    size_t step_size = num_high_risk ? num_agents_to_create / num_high_risk
      : num_agents_to_create;
    for (unsigned i = 0; i < num_agents_to_create; ++i) {
      HIVAgent *a = create_hiv_agent(s.context);
      a->sex = ( (i + s.current_iteration) % 2 == 0) ? sim::MALE : sim::FEMALE;
      if( (i + s.current_iteration) % step_size == 0)
	a->risk = 1;
      else
	a->risk = 0;
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
    double decline_rate = 1.0 - s.context("DECLINE_RATE", stage);
    std::string parm = "_PARTIAL_DECLINE_" + std::to_string(stage);
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
	 ++it) {
      HIVAgent *agent = (HIVAgent *) *it;
      if ( agent->hiv == stage) {
	agent->die(s, "DETERMINISTIC");
	++i;
      }
    }
  }
  s.remove_dead_agents();
}


/*
 * E3
 *
 * People move from HIV stage 1 through to HIV infection WHO stage 4,
 *
 */

void diseaseProgressionDeterministic(sim::Simulation &s)
{
  for (int stage = 1; stage < 4; ++stage) {
    unsigned num_agents_to_move = 0;
    double change_rate = 1.0 - s.context("DISEASE_CHANGE_RATE", stage - 1);
    std::string parm = "_PARTIAL_CHANGE_" + std::to_string(stage);
    s.context.set_if_not_set(parm.c_str(), { 0.0 });
    double &partial_change = s.context.get(parm.c_str())[0];
    double change;
    parm = "_AGENTS_" + std::to_string(stage);
    size_t num_agents = s.context(parm.c_str());

    change = std::max(0.0, change_rate * (num_agents + partial_change));
    partial_change = change - std::floor(change);

    if (change < num_agents)
      num_agents_to_move = num_agents - std::floor(change);

    unsigned i = 0;
    for (auto it = s.agents.begin();
	 it != s.agents.end() && i < num_agents_to_move;
	 ++it) {
      HIVAgent *agent = (HIVAgent *) *it;
      if ( agent->hiv == stage) {
	++agent->hiv;
	++i;
      }
    }
  }
}

/*
 * E4
 *
 * Two risk groups modelled, high and low, adapted (and simplified) from pages
 * 260-261 of Vynnycky & White: Infectious Disease Modelling, Oxford University
 * Press, 2010.
 *
 */

void infectionRiskDeterministic(sim::Simulation &s)
{
  size_t total_infections = 0;
  s.context.set_if_not_set("_PARTIAL_INFECTIONS_HIGH", {0.0});
  double &partial_infections_high = s.context.get("_PARTIAL_INFECTIONS_HIGH")[0];
  s.context.set_if_not_set("_PARTIAL_INFECTIONS_LOW", {0.0});
  double &partial_infections_low = s.context.get("_PARTIAL_INFECTIONS_LOW")[0];
  size_t infectious = s.context("_AGENTS_1") + s.context("_AGENTS_2") +
    s.context("_AGENTS_3") + s.context("_AGENTS_4");
  double risk_infectious_partner = infectious / s.context("_AGENTS");
  double time_period = s.context("INFECTION_PARAMETERS_TIME");
  double force_infection = s.context("FORCE_INFECTION");
  double risk_individual_act = risk_infectious_partner * force_infection;
  double risk_not_infected_individual_act = 1.0 - risk_individual_act;
  // High risk
  double num_partnerships_high = s.context("NUM_PARTNERSHIPS_HIGH");
  double risk_not_infected_high = pow(risk_not_infected_individual_act,
				      num_partnerships_high);
  double risk_infected_high = 1.0 - risk_not_infected_high;
  risk_infected_high = sim::time_correct_prob(risk_infected_high,
					      time_period,
					      s.time_step);
  double infections_high = (1.0 + risk_infected_high) *
    (s.context("_AGENTS_0_HIGH") + partial_infections_high) -
    s.context("_AGENTS_0_HIGH");
  // Low risk
  double num_partnerships_low = s.context("NUM_PARTNERSHIPS_LOW");
  double risk_not_infected_low = pow(risk_not_infected_individual_act,
				     num_partnerships_low);
  double risk_infected_low = 1.0 - risk_not_infected_low;
  risk_infected_low = sim::time_correct_prob(risk_infected_low,
					     time_period,
					     s.time_step);
  double infections_low = (1.0 + risk_infected_low) *
    (s.context("_AGENTS_0_LOW") + partial_infections_low) -
    s.context("_AGENTS_0_LOW");

  partial_infections_high = infections_high - floor(infections_high);

  partial_infections_low = infections_low - floor(infections_low);

  size_t num_infections_high = floor(infections_high);
  size_t num_infections_low = floor(infections_low);

  for (auto a = s.agents.begin(); a != s.agents.end() && num_infections_high;
       ++a) {
    HIVAgent *agent = (HIVAgent *) *a;
    if (agent->risk == 1 && agent->hiv == 0) {
      --num_infections_high;
      agent->hiv = 1;
      agent->hiv_infection_date = s.current_date;
      ++total_infections;
    }
  }

  for (auto a = s.agents.begin(); a != s.agents.end() && num_infections_low;
       ++a) {
    HIVAgent *agent = (HIVAgent *) *a;
    if (agent->risk == 0 && agent->hiv == 0) {
      --num_infections_low;
      agent->hiv = 1;
      agent->hiv_infection_date = s.current_date;
      ++total_infections;
    }
  }
  s.context.tracking.new_infections.push_back(total_infections);
}

/*
 * E5
 *
 * People in each HIV stage move to ARVs with differing probabilities.  People
 * on ARVs can also stop treatment or have viral rebound and move back to stage
 * 3.
 *
 * ARV_RATE 0, 1, 2 and 3 are the rates of moving onto ARVs for WHO stages 1, 2,
 * 3 and 4 respectively. ARV_RATE 5 is the rate for moving off ARVs.
 */

void arvsDeterministic(sim::Simulation &s)
{
  for (int stage = 1; stage < 6; ++stage) {
    unsigned num_agents_to_move = 0;
    double change_rate = 1.0 - s.context("ARV_RATE", stage - 1);
    std::string parm = "_PARTIAL_ARVS_" + std::to_string(stage);
    s.context.set_if_not_set(parm.c_str(), { 0.0 });
    double &partial_change = s.context.get(parm.c_str())[0];
    double change;
    parm = "_AGENTS_" + std::to_string(stage);
    size_t num_agents = s.context(parm.c_str());

    change = std::max(0.0, change_rate * (num_agents + partial_change));
    partial_change = change - std::floor(change);

    if (change < num_agents)
      num_agents_to_move = num_agents - std::floor(change);

    unsigned i = 0;
    for (auto it = s.agents.begin();
	 it != s.agents.end() && i < num_agents_to_move;
	 ++it) {
      HIVAgent *agent = (HIVAgent *) *it;
      if ( agent->hiv == stage) {
	if (stage < 5)
	  agent->hiv = 5;
	else
	  agent->hiv = 3;
	++i;
      }
    }
  }
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

void createDeterministicAgents(Simulation & s)
{
  double num_agents = s.context("NUM_AGENTS");
  // Prevalence by stage
  double prevalence[] = {s.context("HIV_PREVALENCE", 0),
			 s.context("HIV_PREVALENCE", 1),
			 s.context("HIV_PREVALENCE", 2),
			 s.context("HIV_PREVALENCE", 3)};
  // On ARVs
  double prevalence_arvs = s.context("HIV_PREVALENCE_ARVS");
  double high_risk = s.context("HIGH_RISK_PROPORTION");

  for (size_t i = 0; i < num_agents; ++i) {
    sim::HIVAgent *agent = new sim::HIVAgent(s.context);
    agent->hiv = 0;
    agent->risk = 0;
    agent->sex = (i % 2 == 0) ? MALE : FEMALE;
    s.agents.push_back(agent);
  }

  size_t i = num_agents = 0;
  for (size_t stage = 0; stage < 4; ++stage) {
    num_agents += prevalence[stage] * s.agents.size();
    for (; i < num_agents && i < s.agents.size(); ++i)
      ( (HIVAgent *) s.agents[i])->hiv = stage + 1;
  }

  num_agents += prevalence_arvs * s.agents.size();
  for (; i < num_agents && i < s.agents.size(); ++i)
    ( (HIVAgent *) s.agents[i])->hiv = 5;


  num_agents = high_risk * s.agents.size();

  for (i = 0; i < num_agents && i < s.agents.size(); ++i)
    ( (HIVAgent *) s.agents[i])->risk = 1;

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
  unsigned num_stage_1 = count_if(s.agents.begin(), s.agents.end(),
			       [](const sim::Agent *a) {
				 const HIVAgent *b = (HIVAgent *) a;
				 return b->hiv == 1;
			       });
  unsigned num_stage_2 = count_if(s.agents.begin(), s.agents.end(),
			       [](const sim::Agent *a) {
				 const HIVAgent *b = (HIVAgent *) a;
				 return b->hiv == 2;
			       });
  unsigned num_stage_3 = count_if(s.agents.begin(), s.agents.end(),
			       [](const sim::Agent *a) {
				 const HIVAgent *b = (HIVAgent *) a;
				 return b->hiv == 3;
			       });
  unsigned num_stage_4 = count_if(s.agents.begin(), s.agents.end(),
			       [](const sim::Agent *a) {
				 const HIVAgent *b = (HIVAgent *) a;
				 return b->hiv == 4;
			       });
  unsigned num_arvs = count_if(s.agents.begin(), s.agents.end(),
			       [](const sim::Agent *a) {
				 const HIVAgent *b = (HIVAgent *) a;
				 return b->hiv == 5;
			       });
  unsigned num_low_risk_hiv_neg =
    count_if(s.agents.begin(), s.agents.end(),
	     [](const sim::Agent *a) {
	       const HIVAgent *b = (HIVAgent *) a;
	       return b->hiv == 0 && b->risk == 0;
	     });
  unsigned num_low_risk_hiv_pos =
    count_if(s.agents.begin(), s.agents.end(),
	     [](const sim::Agent *a) {
	       const HIVAgent *b = (HIVAgent *) a;
	       return b->hiv > 0 && b->risk == 0;
	     });
  unsigned num_high_risk_hiv_neg =
    count_if(s.agents.begin(), s.agents.end(),
	     [](const sim::Agent *a) {
	       const HIVAgent *b = (HIVAgent *) a;
	       return b->hiv == 0 && b->risk == 1;
	     });
  unsigned num_high_risk_hiv_pos =
    count_if(s.agents.begin(), s.agents.end(),
	     [](const sim::Agent *a) {
	       const HIVAgent *b = (HIVAgent *) a;
	       return b->hiv > 0 && b->risk == 1;
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
  ss << s.simulation_num << ", " << s.current_date << ", stage_1, "
     << num_stage_1 << std::endl;
  ss << s.simulation_num << ", " << s.current_date << ", stage_2, "
     << num_stage_2 << std::endl;
  ss << s.simulation_num << ", " << s.current_date << ", stage_3, "
     << num_stage_3 << std::endl;
  ss << s.simulation_num << ", " << s.current_date << ", stage_4, "
     << num_stage_4 << std::endl;
  ss << s.simulation_num << ", " << s.current_date << ", low_risk_hiv_neg, "
     << num_low_risk_hiv_neg << std::endl;
  ss << s.simulation_num << ", " << s.current_date << ", low_risk_hiv_pos, "
     << num_low_risk_hiv_pos << std::endl;
  ss << s.simulation_num << ", " << s.current_date << ", high_risk_hiv_neg, "
     << num_high_risk_hiv_neg << std::endl;
  ss << s.simulation_num << ", " << s.current_date << ", high_risk_hiv_pos, "
     << num_high_risk_hiv_pos << std::endl;
  ss << s.simulation_num << ", " << s.current_date << ", moved_out, "
     << num_moved_out << std::endl;

  std::cout << ss.str();
}

void testSti(tst::TestSeries &t)
{
  // Test of E1
  std::cout << "Testing E1 time step 1" << std::endl;
  sim::Simulation(sim::Options()
		  .events({sim::advanceTimeEvent
			, calcVariablesEvent
			, increasePopulationDeterministic})
		  .afterEachSimulation([&t](sim::Simulation &s) {
		      TESTEQ(t, s.agents.size(), 1220, "Population growth");
		    })
		  .agentCreate(create_hiv_agent)
		  .timeAdjust("GROWTH", 1.0, 0, 1, sim::PROBABILITY)
		  .parameter("NUM_AGENTS", {1000.0} )
		  .parameter("HIGH_RISK_PROPORTION", {0.5 } )
		  .parameter("TIME_STEP", {1.0} )
		  .parameter("INITIAL_AGE", {15.0} )
		  .parameter("GROWTH", {0.01} )
		  .parameter("RISK_SUSCEPTIBLE", {0.0} )
		  .parameter("RISK_SUSCEPTIBLE_STDEV", {0.0} )
		  .parameter("RISK_UNINFECTIOUS", {0.0} )
		  .parameter("RISK_UNINFECTIOUS_STDEV", {0.0} )
		  .parameter("PARTNER_FORMATION_RATE", {0.0} )
		  .parameter("PARTNER_FORMATION_RATE_STDEV", {0.0} )
		  .parameter("HIGH_RISK_PROPORTION", {0.5 } )
		  .parameter("PREVALENCE", {0.0} )).simulate();

  std::cout << "Testing E1 time step half year" << std::endl;
  sim::Simulation(sim::Options()
		  .events({sim::advanceTimeEvent
			, calcVariablesEvent
			, increasePopulationDeterministic})
		  .afterEachSimulation([&t](sim::Simulation &s) {
		      TESTEQ(t, s.agents.size(), 1221, "Population growth");
		    })
		  .agentCreate(create_hiv_agent)
		  .timeAdjust("GROWTH", 1.0, 0, 1, sim::PROBABILITY)
		  .parameter("NUM_AGENTS", {1000.0} )
		  .parameter("HIGH_RISK_PROPORTION", {0.5 } )
		  .parameter("TIME_STEP", {0.5} )
		  .parameter("INITIAL_AGE", {15.0} )
		  .parameter("GROWTH", {0.01} )
		  .parameter("RISK_SUSCEPTIBLE", {0.0} )
		  .parameter("RISK_SUSCEPTIBLE_STDEV", {0.0} )
		  .parameter("RISK_UNINFECTIOUS", {0.0} )
		  .parameter("RISK_UNINFECTIOUS_STDEV", {0.0} )
		  .parameter("PARTNER_FORMATION_RATE", {0.0} )
		  .parameter("PARTNER_FORMATION_RATE_STDEV", {0.0} )
		  .parameter("HIGH_RISK_PROPORTION", {0.5 } )
		  .parameter("PREVALENCE", {0.0} )).simulate();

  std::cout << "Testing E1 time step 1 week" << std::endl;
  sim::Simulation(sim::Options()
		  .events({sim::advanceTimeEvent
			, calcVariablesEvent
			, increasePopulationDeterministic})
		  .afterEachSimulation([&t](sim::Simulation &s) {
		      TESTEQ(t, s.agents.size(), 1222, "Population growth");
		    })
		  .agentCreate(create_hiv_agent)
		  .timeAdjust("GROWTH", 1.0, 0, 1, sim::PROBABILITY)
		  .parameter("NUM_AGENTS", {1000.0} )
		  .parameter("HIGH_RISK_PROPORTION", {0.5 } )
		  .parameter("TIME_STEP", {sim:WEEK} )
		  .parameter("INITIAL_AGE", {15.0} )
		  .parameter("GROWTH", {0.01} )
		  .parameter("RISK_SUSCEPTIBLE", {0.0} )
		  .parameter("RISK_SUSCEPTIBLE_STDEV", {0.0} )
		  .parameter("RISK_UNINFECTIOUS", {0.0} )
		  .parameter("RISK_UNINFECTIOUS_STDEV", {0.0} )
		  .parameter("PARTNER_FORMATION_RATE", {0.0} )
		  .parameter("PARTNER_FORMATION_RATE_STDEV", {0.0} )
		  .parameter("HIGH_RISK_PROPORTION", {0.5 } )
		  .parameter("PREVALENCE", {0.0} )).simulate();


  // Test E2 time step 1
  std::cout << "Testing E2 time step 1" << std::endl;
  sim::Simulation(sim::Options()
		  .events({sim::advanceTimeEvent
			, calcVariablesEvent
			, decreasePopulationDeterministic})
		  .afterEachSimulation([&t](sim::Simulation &s) {
		      TESTEQ(t, s.agents.size(), 817, "Population decline");
		    })
		  .agentCreate(create_hiv_agent)
		  .timeAdjust("DECLINE_RATE", 1.0, 0, 1, sim::PROBABILITY)
		  .parameter("NUM_AGENTS", {1000.0} )
		  .parameter("HIGH_RISK_PROPORTION", {0.5 } )
		  .parameter("DECLINE_RATE",
			     {0.01, 0.02, 0.03, 0.1, 0.5, 0.015} )
		  .parameter("TIME_STEP", {1.0} )
		  .parameter("INITIAL_AGE", {15.0} )
		  .parameter("HIV_PREVALENCE", {0.0, 0.0, 0.0, 0.0} )
		  .parameter("HIGH_RISK_PROPORTION", {0.5 } )).simulate();

  // Test E2 time step half-year
  std::cout << "Testing E2 time step half year" << std::endl;
  sim::Simulation(sim::Options()
		  .events({sim::advanceTimeEvent
			, calcVariablesEvent
			, decreasePopulationDeterministic})
		  .afterEachSimulation([&t](sim::Simulation &s) {
		      TESTEQ(t, s.agents.size(), 817, "Population decline");
		    })
		  .agentCreate(create_hiv_agent)
		  .timeAdjust("DECLINE_RATE", 1.0, 0, 1, sim::PROBABILITY)
		  .parameter("NUM_AGENTS", {1000.0} )
		  .parameter("HIGH_RISK_PROPORTION", {0.5 } )
		  .parameter("DECLINE_RATE",
			     {0.01, 0.02, 0.03, 0.1, 0.5, 0.015} )
		  .parameter("TIME_STEP", { 0.5 } )
		  .parameter("INITIAL_AGE", {15.0} )
		  .parameter("HIV_PREVALENCE", {0.0, 0.0, 0.0, 0.0} )
		  .parameter("HIGH_RISK_PROPORTION", {0.5 } )).simulate();

  // Test E2 time step week
  std::cout << "Testing E2 time step week" << std::endl;
  sim::Simulation(sim::Options()
		  .events({sim::advanceTimeEvent
			, calcVariablesEvent
			, decreasePopulationDeterministic})
		  .afterEachSimulation([&t](sim::Simulation &s) {
		      TESTEQ(t, s.agents.size(), 817, "Population decline");
		    })
		  .agentCreate(create_hiv_agent)
		  .timeAdjust("DECLINE_RATE", 1.0, 0, 1, sim::PROBABILITY)
		  .parameter("NUM_AGENTS", {1000.0} )
		  .parameter("HIGH_RISK_PROPORTION", {0.5 } )
		  .parameter("DECLINE_RATE",
			     {0.01, 0.02, 0.03, 0.1, 0.5, 0.015} )
		  .parameter("TIME_STEP", { sim::WEEK } )
		  .parameter("INITIAL_AGE", {15.0} )
		  .parameter("HIV_PREVALENCE", {0.0, 0.0, 0.0, 0.0} )
		  .parameter("HIGH_RISK_PROPORTION", {0.5 } )).simulate();


  // Test E3
  std::cout << "Testing E3 time step year" << std::endl;
  sim::Simulation(sim::Options()
		  .events({sim::advanceTimeEvent
			, calcVariablesEvent
			, diseaseProgressionDeterministic})
		  .afterEachSimulation([&t](sim::Simulation &s) {
		      unsigned stage_1 =
			count_if(s.agents.begin(), s.agents.end(),
				 [](const sim::Agent *a) {
				   const HIVAgent *b = (HIVAgent *) a;
				   return b->hiv == 1;
				 });
		      TESTEQ(t, stage_1, 11, "Infection change");
		    })
		  .agentCreate(create_hiv_agent)
		  .timeAdjust("DISEASE_CHANGE_RATE", 1.0, 0, 1,
			      sim::PROBABILITY)
		  .parameter("NUM_AGENTS", {1000.0} )
		  .parameter("HIGH_RISK_PROPORTION", {0.5 } )
		  .parameter("DISEASE_CHANGE_RATE", {0.2, 0.2, 0.2, 0.3} )
		  .parameter("TIME_STEP", {1.0 } )
		  .parameter("INITIAL_AGE", {15.0} )
		  .parameter("HIV_PREVALENCE", {1.0, 0.0, 0.0, 0.0} )
		  .parameter("HIGH_RISK_PROPORTION", {0.5 } )).simulate();

  std::cout << "Testing E3 time step half-year" << std::endl;
  sim::Simulation(sim::Options()
		  .events({sim::advanceTimeEvent
			, calcVariablesEvent
			, diseaseProgressionDeterministic})
		  .afterEachSimulation([&t](sim::Simulation &s) {
		      unsigned stage_1 =
			count_if(s.agents.begin(), s.agents.end(),
				 [](const sim::Agent *a) {
				   const HIVAgent *b = (HIVAgent *) a;
				   return b->hiv == 1;
				 });
		      TESTEQ(t, stage_1, 11, "Infection change");
		    })
		  .agentCreate(create_hiv_agent)
		  .timeAdjust("DISEASE_CHANGE_RATE", 1.0, 0, 1,
			      sim::PROBABILITY)
		  .parameter("NUM_AGENTS", {1000.0} )
		  .parameter("HIGH_RISK_PROPORTION", {0.5 } )
		  .parameter("DISEASE_CHANGE_RATE", {0.2, 0.2, 0.2, 0.3} )
		  .parameter("TIME_STEP", {0.5 } )
		  .parameter("INITIAL_AGE", {15.0} )
		  .parameter("HIV_PREVALENCE", {1.0, 0.0, 0.0, 0.0} )
		  .parameter("HIGH_RISK_PROPORTION", {0.5 } )).simulate();

  std::cout << "Testing E3 time step week." << std::endl;
  sim::Simulation(sim::Options()
		  .events({sim::advanceTimeEvent
			, calcVariablesEvent
			, diseaseProgressionDeterministic})
		  .afterEachSimulation([&t](sim::Simulation &s) {
		      unsigned stage_1 =
			count_if(s.agents.begin(), s.agents.end(),
				 [](const sim::Agent *a) {
				   const HIVAgent *b = (HIVAgent *) a;
				   return b->hiv == 1;
				 });
		      TESTEQ(t, stage_1, 11, "Infection change");
		    })
		  .agentCreate(create_hiv_agent)
		  .timeAdjust("DISEASE_CHANGE_RATE", 1.0, 0, 1,
			      sim::PROBABILITY)
		  .parameter("NUM_AGENTS", {1000.0} )
		  .parameter("HIGH_RISK_PROPORTION", {0.5 } )
		  .parameter("DISEASE_CHANGE_RATE", {0.2, 0.2, 0.2, 0.3} )
		  .parameter("TIME_STEP", {sim::WEEK} )
		  .parameter("INITIAL_AGE", {15.0} )
		  .parameter("HIV_PREVALENCE", {1.0, 0.0, 0.0, 0.0} )
		  .parameter("HIGH_RISK_PROPORTION", {0.5 } )).simulate();



  // Test E4 - I can't figure out how to time correct these parameters. The
  // equations are too hard. See calculations.ods.

  std::cout << "Testing E4 one year with none on ARVs" << std::endl;
  sim::Simulation(sim::Options()
		  .events({sim::advanceTimeEvent
			, calcVariablesEvent
			, infectionRiskDeterministic})
		  .afterEachSimulation([&t](sim::Simulation &s) {
		      unsigned susceptible =
			count_if(s.agents.begin(), s.agents.end(),
				 [](const sim::Agent *a) {
				   const HIVAgent *b = (HIVAgent *) a;
				   return b->hiv == 0;
				 });
		      TESTEQ(t, susceptible, 589, "Infections");
		    })
		  .allAgentsCreate(createDeterministicAgents)
		  .parameter("NUM_AGENTS", {1000.0} )
		  .parameter("HIGH_RISK_PROPORTION", {0.5 } )
		  .parameter("TIME_STEP", {1.0} )
		  .parameter("INITIAL_AGE", {15.0} )
		  .parameter("HIGH_RISK_PROPORTION", {0.5 } )
		  .parameter("INFECTION_PARAMETERS_TIME", {1.0})
		  .parameter("FORCE_INFECTION", {0.05} )
		  .parameter("NUM_PARTNERSHIPS_LOW", {0} )
		  .parameter("NUM_PARTNERSHIPS_HIGH", {6} )
		  .parameter("HIV_PREVALENCE", {0.1, 0.00, 0.00, 0.00} )
		  .parameter("HIV_PREVALENCE_ARVS", {0.00} ))
		  .simulate();

  std::cout << "Testing E4 one year with some on ARVs" << std::endl;
  sim::Simulation(sim::Options()
		  .events({sim::advanceTimeEvent
			, calcVariablesEvent
			, infectionRiskDeterministic})
		  .afterEachSimulation([&t](sim::Simulation &s) {
		      unsigned susceptible =
			count_if(s.agents.begin(), s.agents.end(),
				 [](const sim::Agent *a) {
				   const HIVAgent *b = (HIVAgent *) a;
				   return b->hiv == 0;
				 });
		      TESTEQ(t, susceptible, 664, "Infections");
		    })
		  .allAgentsCreate(createDeterministicAgents)
		  .parameter("NUM_AGENTS", {1000.0} )
		  .parameter("HIGH_RISK_PROPORTION", {0.5 } )
		  .parameter("TIME_STEP", {1.0} )
		  .parameter("INITIAL_AGE", {15.0} )
		  .parameter("HIGH_RISK_PROPORTION", {0.5 } )
		  .parameter("INFECTION_PARAMETERS_TIME", {1.0})
		  .parameter("FORCE_INFECTION", {0.05} )
		  .parameter("NUM_PARTNERSHIPS_LOW", {0} )
		  .parameter("NUM_PARTNERSHIPS_HIGH", {6} )
		  .parameter("HIV_PREVALENCE", {0.05, 0.00, 0.00, 0.00} )
		  .parameter("HIV_PREVALENCE_ARVS", {0.1} ))
		  .simulate();

  std::cout << "Testing E4 one year with a lot on ARVs" << std::endl;
  sim::Simulation(sim::Options()
		  .events({sim::advanceTimeEvent
			, calcVariablesEvent
			, infectionRiskDeterministic})
		  .afterEachSimulation([&t](sim::Simulation &s) {
		      unsigned susceptible =
			count_if(s.agents.begin(), s.agents.end(),
				 [](const sim::Agent *a) {
				   const HIVAgent *b = (HIVAgent *) a;
				   return b->hiv == 0;
				 });
		      TESTEQ(t, susceptible, 576, "Infections");
		    })
		  .allAgentsCreate(createDeterministicAgents)
		  .parameter("NUM_AGENTS", {1000.0} )
		  .parameter("HIGH_RISK_PROPORTION", {0.5 } )
		  .parameter("TIME_STEP", {1.0} )
		  .parameter("INITIAL_AGE", {15.0} )
		  .parameter("HIGH_RISK_PROPORTION", {0.5 } )
		  .parameter("INFECTION_PARAMETERS_TIME", {1.0})
		  .parameter("FORCE_INFECTION", {0.05} )
		  .parameter("NUM_PARTNERSHIPS_LOW", {0} )
		  .parameter("NUM_PARTNERSHIPS_HIGH", {6} )
		  .parameter("HIV_PREVALENCE", {0.1, 0.00, 0.00, 0.00} )
		  .parameter("HIV_PREVALENCE_ARVS", {0.2} ))
		  .simulate();



  // Test E5 - still need to verify this on a spreadsheet
  std::cout << "Testing E5 one year - ARVS to WHO stage only" << std::endl;
  sim::Simulation(sim::Options()
  		  .events({sim::advanceTimeEvent
  			, calcVariablesEvent
  			, arvsDeterministic})
  		  .afterEachSimulation([&t](sim::Simulation &s) {
  		      unsigned on_arvs =
  			count_if(s.agents.begin(), s.agents.end(),
  				 [](const sim::Agent *a) {
  				   const HIVAgent *b = (HIVAgent *) a;
  				   return b->hiv == 5;
  				 });
  		      TESTEQ(t, on_arvs, 60, "On arvs");
  		    })
  		  .allAgentsCreate(createDeterministicAgents)
  		  .timeAdjust("ARV_RATE", 1.0, 0, 1, sim::PROBABILITY)
  		  .parameter("ARV_RATE", {0.0, 0.0, 0.0, 0.0, 0.1} )
  		  .parameter("NUM_AGENTS", {1000.0} )
  		  .parameter("TIME_STEP", {1.0 } )
  		  .parameter("HIV_PREVALENCE", {0.0, 0.0, 0.0, 0.0} )
  		  .parameter("HIV_PREVALENCE_ARVS", {0.5} )
  		  .parameter("HIGH_RISK_PROPORTION", {0.5 } )
  		  .parameter("INITIAL_AGE", {15.0} )).simulate();

  std::cout << "Testing E5 week - ARVS to WHO stage only" << std::endl;
  sim::Simulation(sim::Options()
  		  .events({sim::advanceTimeEvent
  			, calcVariablesEvent
  			, arvsDeterministic})
  		  .afterEachSimulation([&t](sim::Simulation &s) {
  		      unsigned on_arvs =
  			count_if(s.agents.begin(), s.agents.end(),
  				 [](const sim::Agent *a) {
  				   const HIVAgent *b = (HIVAgent *) a;
  				   return b->hiv == 5;
  				 });
  		      TESTEQ(t, on_arvs, 60, "On arvs");
  		    })
  		  .allAgentsCreate(createDeterministicAgents)
  		  .timeAdjust("ARV_RATE", 1.0, 0, 1, sim::PROBABILITY)
  		  .parameter("ARV_RATE", {0.0, 0.0, 0.0, 0.0, 0.1} )
  		  .parameter("NUM_AGENTS", {1000.0} )
  		  .parameter("TIME_STEP", { sim::WEEK } )
  		  .parameter("HIV_PREVALENCE", {0.0, 0.0, 0.0, 0.0} )
  		  .parameter("HIV_PREVALENCE_ARVS", {0.5} )
  		  .parameter("HIGH_RISK_PROPORTION", {0.5 } )
  		  .parameter("INITIAL_AGE", {15.0} )).simulate();

  std::cout << "Testing E5 year - WHO stage only to ARVs" << std::endl;
  sim::Simulation(sim::Options()
  		  .events({sim::advanceTimeEvent
  			, calcVariablesEvent
  			, arvsDeterministic})
  		  .afterEachSimulation([&t](sim::Simulation &s) {
  		      unsigned on_arvs =
  			count_if(s.agents.begin(), s.agents.end(),
  				 [](const sim::Agent *a) {
  				   const HIVAgent *b = (HIVAgent *) a;
  				   return b->hiv == 5;
  				 });
  		      TESTEQ(t, on_arvs, 786, "On arvs");
  		    })
  		  .allAgentsCreate(createDeterministicAgents)
  		  .timeAdjust("ARV_RATE", 1.0, 0, 1, sim::PROBABILITY)
  		  .parameter("ARV_RATE", {0.1, 0.2, 0.3, 0.5, 0.0} )
  		  .parameter("NUM_AGENTS", {1000.0} )
  		  .parameter("TIME_STEP", { sim::YEAR } )
  		  .parameter("HIV_PREVALENCE", {0.1, 0.2, 0.3, 0.1} )
  		  .parameter("HIV_PREVALENCE_ARVS", {0.1} )
  		  .parameter("HIGH_RISK_PROPORTION", {0.5 } )
  		  .parameter("INITIAL_AGE", {15.0} )).simulate();

  std::cout << "Testing E5 year - WHO stage only to ARVs" << std::endl;
  sim::Simulation(sim::Options()
  		  .events({sim::advanceTimeEvent
  			, calcVariablesEvent
  			, arvsDeterministic})
  		  .afterEachSimulation([&t](sim::Simulation &s) {
  		      unsigned on_arvs =
  			count_if(s.agents.begin(), s.agents.end(),
  				 [](const sim::Agent *a) {
  				   const HIVAgent *b = (HIVAgent *) a;
  				   return b->hiv == 5;
  				 });
  		      TESTEQ(t, on_arvs, 786, "On arvs");
  		    })
  		  .allAgentsCreate(createDeterministicAgents)
  		  .timeAdjust("ARV_RATE", 1.0, 0, 1, sim::PROBABILITY)
  		  .parameter("ARV_RATE", {0.1, 0.2, 0.3, 0.5, 0.0} )
  		  .parameter("NUM_AGENTS", {1000.0} )
  		  .parameter("TIME_STEP", { sim::WEEK } )
  		  .parameter("HIV_PREVALENCE", {0.1, 0.2, 0.3, 0.1} )
  		  .parameter("HIV_PREVALENCE_ARVS", {0.1} )
  		  .parameter("HIGH_RISK_PROPORTION", {0.5 } )
  		  .parameter("INITIAL_AGE", {15.0} )).simulate();
}




int main(int argc, char **argv)
{

  // Model A

  sim::Simulation(sim::Options()
		  .additionalTests({testSti})
		  .allAgentsCreate(createDeterministicAgents)
		  .beforeEachSimulation(report)
		  .events({sim::advanceTimeEvent
			, calcVariablesEvent
			, diseaseProgressionDeterministic
			, arvsDeterministic
			, infectionRiskDeterministic
			, decreasePopulationDeterministic
			, increasePopulationDeterministic
			, report
			})
		  .timeAdjust("GROWTH", 1.0, 0, 1, sim::PROBABILITY)
		  .timeAdjust("DISEASE_CHANGE_RATE", 1.0, 0, 1, sim::PROBABILITY)
		  .parameter("INITIAL_AGE", {15.0} )
		  .parameter("GROWTH", {0.02} )
		  .parameter("HIV_PREVALENCE", {0.02, 0.02, 0.02, 0.02} )
		  .parameter("HIV_PREVALENCE_ARVS", {0.02 } )
		  .parameter("HIGH_RISK_PROPORTION", {0.5 } )
  		  .timeAdjust("ARV_RATE", 1.0, 0, 1, sim::PROBABILITY)
  		  .parameter("ARV_RATE", {0.1, 0.2, 0.5, 0.6, 0.2} )
		  .timeAdjust("DECLINE_RATE", 1.0, 0, 1, sim::PROBABILITY)
		  .parameter("DECLINE_RATE",
			     {0.01, 0.02, 0.03, 0.1, 0.5, 0.015} )
		  .parameter("DISEASE_CHANGE_RATE", {0.01, 0.2, 0.2, 0.2} )
		  .parameter("FORCE_INFECTION", {0.05} )
		  .timeAdjust("NUM_PARTNERSHIPS_LOW", 1.0, 0, 1, sim::LINEAR)
		  .parameter("NUM_PARTNERSHIPS_LOW", {1} )
		  .timeAdjust("NUM_PARTNERSHIPS_HIGH", 1.0, 0, 1, sim::LINEAR)
		  .parameter("NUM_PARTNERSHIPS_HIGH", {5} )
		  .afterEachSimulation(report)
		  .commandLine(argc, argv)).simulate();
  return 0;
}
