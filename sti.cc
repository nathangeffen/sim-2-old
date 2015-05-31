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
  double beta = s.context("BETA");
  double stdev = s.context("BETA_STDEV");
  std::normal_distribution<double> normal(beta, stdev);
  double initial_age = s.context("INITIAL_AGE");
  s.context.set_if_not_set("NEW_BORNS", { 0.0 });
  double &context_new_borns = s.context.get("NEW_BORNS")[0];
  double new_borns;

  beta = normal(sim::rng);
  new_borns = std::max(0.0, beta * (context_new_borns + s.agents.size()) -
		       s.agents.size());
  num_agents_to_create = (unsigned) std::floor(new_borns);
  context_new_borns = new_borns - floor(new_borns);

  for (unsigned i = 0; i < num_agents_to_create; ++i) {
    StiAgent *a = new StiAgent(s.context);
    a->sex = uni(sim::rng) < s.context("PROB_MALE") ?
			     sim::MALE : sim::FEMALE;
    a->dob = s.current_date - initial_age;
    a->status = SUSCEPTIBLE;
    s.agents.push_back(a);
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

int main(int argc, char **argv)
{
  sim::Simulation(sim::Options()
		  .events({sim::advanceTimeEvent
			, increasePopulation})
		  .afterEachSimulation(report)
		  .commandLine(argc, argv)
		  .agentCreate(createStiAgent)).simulate();
  return 0;
}
