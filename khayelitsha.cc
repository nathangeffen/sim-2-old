#include "sim.hh"

class KhayaAgent : public sim::HIVAgent {
public:
  KhayaAgent(sim::Context &c) : sim::HIVAgent(c) {};
  std::string area;
};

void report(sim::Simulation &s)
{
  unsigned num_alive = s.agents.size();
  unsigned num_dead = s.dead_agents.size();
  unsigned num_moved_out = s.aged_out_agents.size();
  unsigned num_site_b = std::count_if(s.agents.begin(), s.agents.end(),
				      [](Agent *a) {
					KhayaAgent *agent = (KhayaAgent *) a;
					return (agent->area == "Site B");
				      });
  unsigned num_site_c = std::count_if(s.agents.begin(), s.agents.end(),
				      [](Agent *a) {
					KhayaAgent *agent = (KhayaAgent *) a;
					return (agent->area == "Site C");
				      });
  unsigned num_males = std::count_if(s.agents.begin(), s.agents.end(),
				     [](Agent *a) {
				       return (a->sex == sim::MALE);
				     });
  unsigned num_females = std::count_if(s.agents.begin(), s.agents.end(),
				       [](Agent *a) {
					 return (a->sex == sim::FEMALE);
				       });
  unsigned num_under_5 = std::count_if(s.agents.begin(), s.agents.end(),
				       [&s](Agent *a) {
					 return
					 (s.context("START_DATE") - a->dob < 5);
				       });
  unsigned num_5_10 = std::count_if(s.agents.begin(), s.agents.end(),
				    [&s](Agent *a) {
				      return
				      (s.context("START_DATE") - a->dob >= 5) &&
				      (s.context("START_DATE") - a->dob < 10);
				    });

  std::stringstream ss;
  ss << s.simulation_num << ", " << s.current_date << ", alive, "
     << num_alive << std::endl;
  ss << s.simulation_num << ", " << s.current_date << ", dead, "
     << num_dead << std::endl;
  ss << s.simulation_num << ", " << s.current_date << ", moved_out, "
     << num_moved_out << std::endl;
  ss << s.simulation_num << ", " << s.current_date << ", Site B, "
     << num_site_b << std::endl;
  ss << s.simulation_num << ", " << s.current_date << ", Site C, "
     << num_site_c << std::endl;
  ss << s.simulation_num << ", " << s.current_date << ", Males, "
     << num_males << std::endl;
  ss << s.simulation_num << ", " << s.current_date << ", Females, "
     << num_females << std::endl;
  ss << s.simulation_num << ", " << s.current_date << ", Under 5, "
     << num_under_5 << std::endl;
  ss << s.simulation_num << ", " << s.current_date << ", 5-10, "
     << num_5_10 << std::endl;
  std::cout << ss.str();
}

int numFrom(std::string val)
{
  return stoi(val.substr(0, val.find("-")));
}

int numTo(std::string val)
{
  return stoi(val.substr(val.find("-") + 1));
}

void createAgents(Simulation & s)
{
  double start_date = s.context("START_DATE");
  for (auto &i : s.context.get_initial()) {
    for (size_t j = 0; j < i.first; j++) {
      auto & fields = i.second;
      KhayaAgent *a = new KhayaAgent(s.context);
      {
	int age_lower = numFrom(fields[0]);
	int age_higher = numTo(fields[0]);

	std::uniform_int_distribution<> dis(age_lower, age_higher);
	a->dob = start_date - dis(sim::rng);
      }
      a->area = fields[1];
      if (fields[2] == "male")
	a->sex = MALE;
      else
	a->sex = FEMALE;
      if (fields[3] == "pos")
	a->hiv = 1;
      else
	a->hiv = 0;
      s.agents.push_back(a);
    }
  }
}

void calcVariablesEvent(sim::Simulation &s)
{
  size_t num_neg = 0.0, num_pos = 0.0, num_arvs = 0.0;

  for (auto & a: s.agents) {
    KhayaAgent *agent = (KhayaAgent *) a;
    if (agent->hiv == 0)
      ++num_neg;
    else {
      ++num_pos;
      if (agent->on_arvs_date > 0.0)
	++num_arvs;
    }
  }
  s.context.set("_AGENTS",s.agents.size());
  s.context.set("_NEG", num_neg);
  s.context.set("_POS", num_pos);
  s.context.set("_ARVS", num_arvs);
}


int main(int argc, char **argv)
{
  sim::Simulation(sim::Options()
		  .allAgentsCreate(createAgents)
		  .events({sim::advanceTimeEvent
			, calcVariablesEvent
			// , increasePopulationDiffEqEvent
			// , emigrationDiffEqEvent
			// , infectedDiffEqEvent
			})
		  .afterEachSimulation(report)
		  .commandLine(argc, argv)).simulate();
  return 0;
}
