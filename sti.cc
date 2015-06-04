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

void increasePopulation(sim::Simulation &s)
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

void emigration(sim::Simulation &s)
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

  sim::Simulation(sim::Options()
		  .events({sim::advanceTimeEvent, increasePopulation})
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
  		  .events({sim::advanceTimeEvent, emigration})
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
}

int main(int argc, char **argv)
{
  sim::Simulation(sim::Options()
		  .additionalTests({testSti})
		  .events({sim::advanceTimeEvent //, increasePopulation
			})
		  .afterEachSimulation(report)
		  .commandLine(argc, argv)
		  .agentCreate(createStiAgent)).simulate();
  return 0;
}
