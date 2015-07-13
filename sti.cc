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

enum Status {
  SUSCEPTIBLE = 0,
  INFECTED = 1,
  UNINFECTIOUS = 2
};

class StiAgent : public sim::Agent {
public:
  StiAgent(sim::Context &c) : sim::Agent(c) {};
  Status status;
  std::vector<double> dates_infected;
  double date_last_child = 0.0;
  double moved_out = 0.0;
  double risk_susceptible;
  double risk_uninfectious;
  double rate_partner_change;
};

StiAgent *createStiAgent(sim::Context &c)
{
  std::uniform_real_distribution<double> uni;
  StiAgent *a = new StiAgent(c);

  a->sex = uni(sim::rng) < c("PROB_MALE") ?
			   sim::MALE : sim::FEMALE;
  { // dob
    std::uniform_real_distribution<double>
      uni_age(c("EARLIEST_BIRTH_DATE"),
	      c("LATEST_BIRTH_DATE"));
    a->dob = uni_age(sim::rng);
  }

  { // risk of infection with partner with STI, if susceptible
    double risk_susceptible = c("RISK_SUSCEPTIBLE");
    double risk_susceptible_stdev = c("RISK_SUSCEPTIBLE_STDEV");
    std::normal_distribution<double> risk(risk_susceptible,
					  risk_susceptible_stdev);
    a->risk_susceptible = risk(sim::rng);
  }

  { // risk of infection with partner with STI, if uninfectious
    double risk_uninfectious = c("RISK_UNINFECTIOUS");
    double risk_uninfectious_stdev = c("RISK_UNINFECTIOUS_STDEV");
    std::normal_distribution<double> risk(risk_uninfectious,
					  risk_uninfectious_stdev);
    a->risk_uninfectious = risk(sim::rng);
  }

  { // Rate of partner change
    double rate_partner_change = c("PARTNER_FORMATION_RATE");
    double rate_partner_change_stdev = c("PARTNER_FORMATION_RATE_STDEV");
    std::normal_distribution<double> rate(rate_partner_change,
					  rate_partner_change_stdev);
    a->rate_partner_change = rate(sim::rng);
  }

  // Status
  if (uni(sim::rng) < c("PREVALENCE"))
    a->status = INFECTED;
  else
    a->status = SUSCEPTIBLE;

  return a;
}

void calcVariablesEvent(sim::Simulation &s)
{
  size_t num_susceptible = 0.0, num_infected = 0.0, num_uninfectious = 0.0;

  for (auto & a: s.agents) {
    StiAgent *agent = (StiAgent *) a;
    if (agent->status == SUSCEPTIBLE)
      ++num_susceptible;
    else if (agent->status == INFECTED)
      ++num_infected;
    else
      ++num_uninfectious;

  }
  s.context.set("_AGENTS",s.agents.size());
  s.context.set("_SUSCEPTIBLE", num_susceptible);
  s.context.set("_INFECTIOUS", num_infected);
  s.context.set("_UNINFECTIOUS", num_uninfectious);
}

/* E1 */

void increasePopulationDeterministic(sim::Simulation &s)
{
  std::uniform_real_distribution<double> uni;
  unsigned num_agents_to_create;
  double growth_rate = s.context("GROWTH");
  // double stdev = s.context("GROWTH_STDEV");
  // std::normal_distribution<double> normal(growth_rate, stdev);
  double initial_age = s.context("INITIAL_AGE");
  s.context.set_if_not_set("_PARTIAL_GROWTH", { 0.0 });
  double &partial_growth = s.context.get("_PARTIAL_GROWTH")[0];
  double growth;
  size_t num_agents = s.context("_AGENTS");

  // growth_rate = normal(sim::rng);
  growth = std::max(0.0, growth_rate * (num_agents + partial_growth));
  partial_growth = growth - std::floor(growth);

  if (growth > num_agents) {
    num_agents_to_create = std::floor(growth) - num_agents;

    for (unsigned i = 0; i < num_agents_to_create; ++i) {
      StiAgent *a = createStiAgent(s.context);
      a->sex = (i % 2 == 0) ? sim::MALE : sim::FEMALE;
      a->dob = s.current_date - initial_age;
      a->status = SUSCEPTIBLE;
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
    unsigned num_agents_to_remove;
    std::string parm = "DECLINE_RATE_" + stage;
    double decline_rate = s.context(parm.c_str());
    parm = "_PARTIAL_DECLINE_" + stage;
    s.context.set_if_not_set(parm.c_str(), { 0.0 });
    double &partial_decline = s.context.get(parm.c_str())[0];
    double decline, floor_decline;
    parm = "_AGENTS_" + stage;
    size_t num_agents = s.context(parm.c_str());

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

void increasePopulationStochasticEvent(sim::Simulation &s)
{
  double growth_rate = s.context("GROWTH_STOCHASTIC");
  double initial_age = s.context("INITIAL_AGE");
  size_t num_agents = s.context("_AGENTS");

  for (size_t i = 0; i < num_agents; ++i) {
    auto sti_agent = (StiAgent *) s.agents[i];
    if (s.current_date - sti_agent->date_last_child > 1.0) {
      if (sim::is_event(growth_rate)) {
	sti_agent->date_last_child = s.current_date;
	StiAgent *a = createStiAgent(s.context);
	a->sex = sim::is_event(s.context("PROB_MALE")) ?
	  sim::MALE : sim::FEMALE;
	a->dob = s.current_date - initial_age;
	a->status = SUSCEPTIBLE;
	s.agents.push_back(a);
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

void emigrationStochasticEvent(sim::Simulation &s)
{
  double decline_rate = s.context("DECLINE_STOCHASTIC");

  for (auto & agent : s.agents) {
    auto sti_agent = (StiAgent *) agent;
    if (sim::is_event(decline_rate)) {
      sti_agent->moved_out = s.current_date;
    }
  }
  for (size_t i = 0; i < s.agents.size(); ++i) {
    auto sti_agent = (StiAgent *) s.agents[i];
    if (sti_agent->moved_out)
      s.move_agent(s.agents, s.aged_out_agents, i);
  }
}

void infectedDiffEqEvent(sim::Simulation &s)
{
  double_t I = s.context("_INFECTIOUS"); // num infectious
  double S = s.context("_SUSCEPTIBLE");
  double N = s.context("_AGENTS"); // num agents
  double recovery_rate = s.time_step / s.context("DURATION_INFECTION");
  double beta = s.context("RISK_SUSCEPTIBLE"); // transmission probability
  double c = s.context("PARTNER_FORMATION_RATE");
  double delta = c * beta * (I / N);

  double estimated_new_infections = delta * S - recovery_rate * I;
  double risk_individual_infected = estimated_new_infections / S;

  // std::cout << "D0: " << delta << " " << I << " " << S << " " << N << " "
  //	    << std::setprecision(6)  << risk_individual_infected << std::endl;


  for (auto &a : s.agents) {
    StiAgent *agent  = (StiAgent *) a;
    if (agent->status == SUSCEPTIBLE) {
      if (sim::is_event(risk_individual_infected)) {
	agent->status = INFECTED;
	agent->dates_infected.push_back(s.current_date);
      }
    }
  }
}



void report(sim::Simulation &s)
{
  unsigned num_alive = s.agents.size();
  unsigned num_susceptible = count_if(s.agents.begin(), s.agents.end(),
				      [](const sim::Agent *a) {
					const StiAgent *b = (StiAgent *) a;
					return b->status == SUSCEPTIBLE;
				      });
  unsigned num_infected = count_if(s.agents.begin(), s.agents.end(),
				   [](const sim::Agent *a) {
				     const StiAgent *b = (StiAgent *) a;
				     return b->status == INFECTED;
				   });
  unsigned num_uninfectious = count_if(s.agents.begin(), s.agents.end(),
				       [](const sim::Agent *a) {
					 const StiAgent *b = (StiAgent *) a;
					 return b->status == UNINFECTIOUS;
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
  ss << s.simulation_num << ", " << s.current_date << ", uninfectious, "
     << num_uninfectious << std::endl;
  ss << s.simulation_num << ", " << s.current_date << ", dead, "
     << num_dead << std::endl;
  ss << s.simulation_num << ", " << s.current_date << ", moved_out, "
     << num_moved_out << std::endl;

  std::cout << ss.str();
}

void testSti(tst::TestSeries &t)
{

  // auto rpt = [&t](sim::Simulation &s) {
  //   std::stringstream ss;
  //   ss << "agents, " << s.simulation_num << ", "
  //   << s.current_date << ", " << s.agents.size() << std::endl;
  //   std::cout << ss.str();
  // };

  sim::Simulation(sim::Options()
		  .events({sim::advanceTimeEvent
			, calcVariablesEvent
			, increasePopulationDeterministic})
		  .afterEachSimulation([&t](sim::Simulation &s) {
		      TESTEQ(t, s.agents.size(), 67274, "Population growth");
		    })
		  .agentCreate(createStiAgent)
		  .timeAdjust("GROWTH", 1.0, 0, 1, sim::COMPOUND)
		  .parameter("NUM_AGENTS", {10000.0} )
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
			, emigrationDiffEqEvent})
		  .afterEachSimulation([&t](sim::Simulation &s) {
		      TESTEQ(t, s.agents.size(), 423, "Population decline");
		    })
		  .agentCreate(createStiAgent)
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
		  .parameter("PREVALENCE", {0.0} )).simulate();

  sim::Simulation(sim::Options()
		  .events({sim::advanceTimeEvent
			, calcVariablesEvent
			, increasePopulationStochasticEvent
			})
		  .afterEachSimulation([&t](sim::Simulation &s) {
		      TEST(t, s.agents.size() > 6116 && s.agents.size() < 7400,
			   "Stochastic growth");
		    })
		  .agentCreate(createStiAgent)
		  .timeAdjust("GROWTH_STOCHASTIC", 1.0, 0, 1, sim::PROBABILITY)
		  .parameter("NUM_AGENTS", {1000.0} )
		  .parameter("TIME_STEP", {0.5} )
		  .parameter("NUM_YEARS", {20.0} )
		  .parameter("INITIAL_AGE", {15.0} )
		  .parameter("GROWTH_STOCHASTIC", {0.1} )
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
			, emigrationStochasticEvent
			})
		  .agentCreate(createStiAgent)
		  .timeAdjust("DECLINE_STOCHASTIC", 1.0, 0, 1, sim::PROBABILITY)
		  .report(sim::numAgents, sim::mean,
			  [&](const double value) {
			    TEST(t, value > 109 && value < 135,
				 "Stochastic decline mean");
			  })
		  .report(sim::numAgents, sim::median,
			  [&](const double value) {
			    TEST(t, value > 109 && value < 135,
				 "Stochastic decline median");
			  })
		  .report(sim::numAgents,
			  [](const std::vector<double> & values) {
			    return *std::min_element(values.cbegin(),
						     values.cend());
			  },
			  [&](const double value) {
			    TEST(t, value < 120, "stochastic decline min");
			  })
		  .report(sim::numAgents,
			  [](const std::vector<double> & values) {
			    return *std::max_element(values.cbegin(),
						     values.cend());
			  },
			  [&](const double value) {
			    TEST(t, value > 130, "stochastic decline max");
			  })
		  .parameter("NUM_AGENTS", {1000.0} )
		  .parameter("TIME_STEP", {0.1} )
		  .parameter("NUM_SIMULATIONS", {10.0} )
		  .parameter("NUM_YEARS", {20.0} )
		  .parameter("INITIAL_AGE", {15.0} )
		  .parameter("DECLINE_STOCHASTIC", {0.1} )
		  .parameter("RISK_SUSCEPTIBLE", {0.0} )
		  .parameter("RISK_SUSCEPTIBLE_STDEV", {0.0} )
		  .parameter("RISK_UNINFECTIOUS", {0.0} )
		  .parameter("RISK_UNINFECTIOUS_STDEV", {0.0} )
		  .parameter("PARTNER_FORMATION_RATE", {0.0} )
		  .parameter("PARTNER_FORMATION_RATE_STDEV", {0.0} )
		  .parameter("PREVALENCE", {0.0} )).simulate();
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
			// , emigrationDiffEqEvent
			// , infectedDiffEqEvent
			// , report
			})
		  .afterEachSimulation(report)
		  .commandLine(argc, argv)
		  .agentCreate(createStiAgent)).simulate();
  return 0;
}
