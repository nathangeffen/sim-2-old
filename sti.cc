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
  double date_last_child = 0.0;
  double moved_out = 0.0;
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

  // Status
  if (uni(sim::rng) < c("PREVALENCE"))
    a->status = INFECTED;
  else
    a->status = SUSCEPTIBLE;

  return a;
}

void increasePopulationDeterministicEvent(sim::Simulation &s)
{
  std::uniform_real_distribution<double> uni;
  unsigned num_agents_to_create;
  double growth_rate = s.context("GROWTH");
  double stdev = s.context("GROWTH_STDEV");
  std::normal_distribution<double> normal(growth_rate, stdev);
  double initial_age = s.context("INITIAL_AGE");
  s.context.set_if_not_set("_PARTIAL_GROWTH", { 0.0 });
  double &partial_growth = s.context.get("_PARTIAL_GROWTH")[0];
  double growth;
  size_t num_agents = s.agents.size();

  growth_rate = normal(sim::rng);
  growth = std::max(0.0, growth_rate * (num_agents + partial_growth));
  partial_growth = growth - std::floor(growth);

  if (growth > num_agents) {
    num_agents_to_create = std::floor(growth) - num_agents;

    for (unsigned i = 0; i < num_agents_to_create; ++i) {
      StiAgent *a = new StiAgent(s.context);
      a->sex = uni(sim::rng) < s.context("PROB_MALE") ?
			       sim::MALE : sim::FEMALE;
      a->dob = s.current_date - initial_age;
      a->status = SUSCEPTIBLE;
      s.agents.push_back(a);
    }
  }
}

void increasePopulationStochasticEvent(sim::Simulation &s)
{
  double growth_rate = s.context("GROWTH_STOCHASTIC");
  double initial_age = s.context("INITIAL_AGE");
  size_t num_agents = s.agents.size();

  for (size_t i = 0; i < num_agents; ++i) {
    auto sti_agent = (StiAgent *) s.agents[i];
    if (s.current_date - sti_agent->date_last_child > 1.0) {
      if (sim::is_event(growth_rate)) {
	sti_agent->date_last_child = s.current_date;
	StiAgent *a = new StiAgent(s.context);
	a->sex = sim::is_event(s.context("PROB_MALE")) ?
	  sim::MALE : sim::FEMALE;
	a->dob = s.current_date - initial_age;
	a->status = SUSCEPTIBLE;
	s.agents.push_back(a);
      }
    }
  }
}



void emigrationDeterministicEvent(sim::Simulation &s)
{
  std::uniform_real_distribution<double> uni;
  unsigned num_agents_to_remove;
  double decline_rate = s.context("DECLINE");
  double stdev = s.context("DECLINE_STDEV");
  std::normal_distribution<double> normal(decline_rate, stdev);
  s.context.set_if_not_set("_PARTIAL_DECLINE", { 0.0 });
  double &partial_decline = s.context.get("_PARTIAL_DECLINE")[0];
  double decline, floor_decline;
  size_t num_agents = s.agents.size();

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

  auto rpt = [&t](sim::Simulation &s) {
    std::stringstream ss;
    ss << "agents, " << s.simulation_num << ", "
    << s.current_date << ", " << s.agents.size() << std::endl;
    std::cout << ss.str();
  };

  sim::Simulation(sim::Options()
		  .events({sim::advanceTimeEvent
			, increasePopulationDeterministicEvent})
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
		  .parameter("PREVALENCE", {0.0} )).simulate();
  sim::Simulation(sim::Options()
  		  .events({sim::advanceTimeEvent
			, emigrationDeterministicEvent})
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
   		  .parameter("PREVALENCE", {0.0} )).simulate();
  sim::Simulation(sim::Options()
		  .events({sim::advanceTimeEvent
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
		  .parameter("PREVALENCE", {0.0} )).simulate();
  sim::Simulation(sim::Options()
		  .events({sim::advanceTimeEvent
			, emigrationStochasticEvent
			})
		  .afterEachSimulation([&t](sim::Simulation &s) {
		      TEST(t, s.agents.size() > 109 && s.agents.size() < 135,
			   "Stochastic decline");
		    })
		  .agentCreate(createStiAgent)
		  .timeAdjust("DECLINE_STOCHASTIC", 1.0, 0, 1, sim::PROBABILITY)
  		  .parameter("NUM_AGENTS", {1000.0} )
		  .parameter("TIME_STEP", {0.1} )
		  .parameter("NUM_SIMULATIONS", {1.0} )
		  .parameter("NUM_YEARS", {20.0} )
		  .parameter("INITIAL_AGE", {15.0} )
		  .parameter("DECLINE_STOCHASTIC", {0.1} )
		  .parameter("PREVALENCE", {0.0} )).simulate();
}


int main(int argc, char **argv)
{
  sim::Simulation(sim::Options()
		  .additionalTests({testSti})
		  .events({sim::advanceTimeEvent
			, increasePopulationDeterministicEvent
			, emigrationDeterministicEvent
			//, report
			})
		  .afterEachSimulation(report)
		  .commandLine(argc, argv)
		  .agentCreate(createStiAgent)).simulate();
  return 0;
}
