#ifndef __SIM_H__
#define __SIM_H__

#include <execinfo.h>
#include <unistd.h>

#include <cassert>
#include <cfloat>
#include <cmath>
#include <cstdio>
#include <cstdlib>


#include <algorithm>
#include <exception>
#include <iomanip>
#include <iostream>
#include <fstream>
#include <functional>
#include <random>
#include <sstream>
#include <string>
#include <thread>
#include <vector>
#include <unordered_map>

#include "test.hh"

namespace sim {

  class Context;
  class Agent;
  class Simulation;

  // Frequently used time periods
  const double YEAR = 1.0;
  const double MONTH = 1.0 / 12.0;
  const double WEEK = 1.0 / 52.0;
  const double DAY = 1.0 / 365.25;

  // Array of event functions types
  typedef std::vector< std::function<void(Simulation &)> > Events;
  // Array of test function types
  typedef std::vector< std::function<void(tst::TestSeries &)> > Tests;

  // The different types of adjustments that can be applied to
  // a parameter.
  enum TimeAdjustMethod {
    PROBABILITY = 0,
    LINEAR = 1,
    COMPOUND = 2
  };

  // Random number generators used throughout the simulation system
  extern thread_local  std::mt19937_64 rng;

  // Generates a random number and returns true if it is less than risk.
  bool
  is_event(const double risk);

  // Makes a linear correction to a parameter of a specified time period to the
  // actual time period of the simulation.
  double
  time_correct_linear(const double parameter_prob,
		      const double parameter_time_period,
		      const double actual_time_period);

  // Makes a correction to a parameter of a specified time period to the
  // actual time period of the simulation, such that:
  // output = parameter ^ (actual_time_period / parameter_time_period).
  double
  time_correct_compound(const double parameter,
			const double parameter_time_period,
			const double actual_time_period);

  // Makes a correction to a probability for a specified time period to the
  // actual time period of the simulation.
  double
  time_correct_prob(const double probability,
		    const double probability_time_period,
		    const double actual_time_period);

  // Used by functions that correct values for a specified time period to the
  // actual time period of the simulation.
  struct TimeAdjust {
    double time_period;
    size_t from;
    size_t step;
    TimeAdjustMethod method;
  };


  // Advances the simulation date by the user defined time simulation time step.
  // Default time step is a week: 365.25/52.0. Almost any useful simulation
  // would include this event.
  void advanceTimeEvent(Simulation &simulation);

  // Default death event - does nothing.
  void deathEvent(Simulation &simulation);

  enum Sex {
    MALE = 0,
    FEMALE = 1
  };

  // Thrown when an invalid parameter is requested.
  class InvalidParameter : public std::exception
  {
  public:
    InvalidParameter(const std::string & s =
		     "Invalid index into parameter array");
  private:
    std::string msg;
    const char * what() const throw();
  };

  // Thrown when invalid user data is encountered.
  class InvalidData : public InvalidParameter
  {
  public:
    InvalidData(const std::string & s = "Invalid user data.") :
      InvalidParameter(s) {}
  };


  class Context {
  public:
    void set_defaults_not_yet_set();
    double operator()(const char * key, size_t i);
    double operator()(const char * c);
    double operator()(const std::string & s);
    double operator()(const std::string &s, size_t i);
    void set(const char *parameter, double value);
    void set(const char *parameter, std::vector<double> values);
    void set(const std::string& parameter, std::vector<double> values);
    void
    set_if_not_set(const char *key, std::vector<double> values);
    std::vector<double> & get(const char *key);
    void push_back(const std::string& parameter, double value);
    size_t size(const char *parameter);
    void
    set_time_adjust(const std::string & parameter,
		    double time_period,
		    size_t from = 0,
		    size_t step = 1,
		    TimeAdjustMethod method = PROBABILITY);
    TimeAdjust get_time_adjust(const char * s);
    double get_time_adjust(std::string & s);
    void set_initial(size_t num_agents,
		     std::vector< std::string > characteristics);
    std::pair<size_t, std::vector<std::string> > & get_initial(size_t i);
    std::vector< std::pair<size_t, std::vector<std::string> > > &
    get_initial();
    void
    convert_probabilities_to_time_period(const char * key,
					 double parameter_time_period,
					 size_t start = 0,
					 size_t step = 1);
    void
    convert_linear_values_to_time_period(const char * key,
					 double parameter_time_period,
					 size_t start = 0,
					 size_t step = 1);
    void
    convert_compound_product_to_time_period(const char * key,
					    double parameter_time_period,
					    size_t start = 0,
					    size_t step = 1);
    void adjust_parameters_to_time_period();
    void print_parameters();
    unsigned last_agent = 0;
    void *user_data = NULL;
  private:
    std::unordered_map<std::string, std::vector<double> > parameters;
    std::unordered_map<std::string,
		       std::tuple<double, size_t, size_t, TimeAdjustMethod> >
    parameters_to_time_adjust;
    std::vector<std::pair<size_t, std::vector<std::string> > > initialization;
  };

  void process_command_line(Context &context,
			    int argc,
			    char **argv,
			    const Tests & tests);
  class Agent
  {
  public:
    Agent(Context &context);
    double age(Simulation &s);
    void die(Simulation &s, const std::string& cause);
    unsigned id;
    bool alive;
    double dob;
    double dod;
    std::string cause;
    Sex sex;
  private:
    Context &context_;
  };

  Agent* create_default_agent(Context & c);

  void threaded_part_of_sim_loop(Simulation & simulation,
				 size_t from, size_t to);

  double numAgents(const Simulation &s);
  double mean(const std::vector<double> & values);
  double median(std::vector<double> & values);

  class Reporter {
  public:
    Reporter(std::function<double(const Simulation &)> getValueFunc,
	     std::function<double(std::vector<double> &)> calcFunc,
	     std::function<void(const double)> outputFunc);
    void setSize(const size_t numVals);
    void getReportValue(const Simulation &s, size_t id);
    void calculate();
    void output();
  private:
    std::function<double(const Simulation &)> getValueFunc_;
    std::function<double(std::vector<double> &)> calcFunc_;
    std::function<void(const double)> outputFunc_;
    std::vector< double > values_;
    double calcVal_;
  };

  class Options {
  public:
    Options();
    Options& events(Events events);
    Options& additionalTests(Tests tests);
    Options& commandLine(int argc, char **argv);
    Options& agentCreate(std::function<Agent *(Context &)>
			 individual_agent_create_func);
    Options& allAgentsCreate(std::function<void (Simulation &)>
			     all_agents_create_func);
    Options& beforeAllSimulations(std::function<void(Simulation &)>
				  before_all_simulations_func);
    Options& beforeEachSimulation(std::function<void(Simulation &)>
				  before_each_simulation_func);
    Options& afterEachSimulation(std::function<void(Simulation &)>
				 after_each_simulation_func);
    Options& parameter(const char * parameter,
		       const std::vector<double> & values);
    Options& timeAdjust(const std::string & parameter,
			double time_period,
			size_t from = 0,
			size_t step = 1,
			TimeAdjustMethod method = PROBABILITY);
    Options& report(std::function<double(const Simulation &)> getValueFunc,
		    std::function<double(std::vector<double> &)> calcFunc,
		    std::function<void(const double)> outputFunc);
  private:
    friend class Simulation;
    Events events_;
    Tests tests_;
    int argc_;
    char **argv_;
    std::function<Agent *(Context &)> individual_agent_create_func_;
    std::function<void(Simulation &)> all_agents_create_func_;
    std::function<void(Simulation &)> before_all_simulations_func_;
    std::function<void(Simulation &)>  before_each_simulation_func_;
    std::function<void(Simulation &)> after_each_simulation_func_;
    Context context_;
    std::vector<Reporter> reporters_;
  };

  class Simulation
  {
  public:
    Simulation();
    Simulation(const Options & options);
    ~Simulation();
    void setOptions(const Options & options);
    void init(unsigned seed = 0);
    void
    init(unsigned seed,
	 Events events,
	 std::function<Agent *(Context &)> agent_create_func,
	 std::function<void (Simulation &)> create_all_agents_func);
    double
    time_correct_prob(const double parameter_prob,
		      const double parameter_time_period);
    void
    move_agent(std::vector <Agent *> &from_agents,
	       std::vector <Agent *> &to_agents,
	       const std::size_t agent_index);
    void remove_dead_agents();
    void shuffle_agents();
    void simulate();
    void simulate_once();
    std::vector<Agent *> agents;
    std::vector<Agent *> aged_out_agents;
    std::vector<Agent *> dead_agents;
    double current_date, time_step;
    unsigned simulation_num = 0, current_iteration, total_iterations;
    Context context;
    void *user_data = NULL;
  private:
    Events events_;
    Tests tests_;
    int argc_;
    char **argv_;
    std::vector<Reporter> *reporters_ = NULL;
    std::function<Agent *(Context &)> agent_create_func_;
    std::function<void (Simulation &)> create_all_agents_func_;
    std::function<void(Simulation &)> before_all_simulations_func_;
    std::function<void(Simulation &)> before_each_simulation_func_;
    std::function<void(Simulation &)> after_each_simulation_func_;
    friend void
    threaded_part_of_sim_loop(Simulation & simulation, size_t from, size_t to)
    {
      for (size_t i = from; i < to; ++i) {
	Simulation s = simulation;
	// The +1 in next line of code is vital. The simulation parent needs to
	// be identified by the destructor to free shared memory, and this is
	// done by checking if simulation_num == 0. If +1 is left out,
	// then two simulations will have simulation_num = 0.
	s.simulation_num = i + 1;
	s.init(i * 23 + 7);
	if (s.before_each_simulation_func_)
	  s.before_each_simulation_func_(s);
	s.simulate_once();
	if (s.after_each_simulation_func_)
	  s.after_each_simulation_func_(s);
	if (s.reporters_)
	  for (auto & r : *(s.reporters_))
	    r.getReportValue(s, i);
      }
    }
  };

  class HIVAgent : public Agent {
  public:
    HIVAgent(Context &c);
    int hiv;
    double hiv_infection_date;
    double on_arvs_date = 0.0;
    double cd4;
    bool circumcised;
    double orientation;
    double riskiness;
    unsigned num_partners;
    std::vector<unsigned> partners;
    friend std::ostream& operator<<(std::ostream& os, const HIVAgent& agent)
    {
      os << " ID: " << agent.id << " Birth: " << agent.dob
	 << " Death: " << agent.dod
	 << " HIV: " << agent.hiv << " HIV infection date: "
	 << agent.hiv_infection_date << " CD4: " << agent.cd4;
      return os;
    }
  };
  HIVAgent *create_hiv_agent(Context &c);
}

////////////////////

// Inline methods

using namespace sim;

// InvalidParameter inline methods

inline InvalidParameter::InvalidParameter(const std::string & s)
{
  msg = s;
}

inline const char*
InvalidParameter::what() const throw()
{
  return msg.c_str();
}

///////////////

// Context inline methods

///////////////

inline double
Context::operator()(const char * c)
{
  return operator()(c, 0);
}

inline double
Context::operator()(const std::string & s) {
  return operator()(s, 0);
}

inline double
Context::operator()(const std::string &s, size_t i)
{
  return operator()(s.c_str(), i);
}

inline void
Context::set(const char *parameter, double value)
{
  parameters[parameter] = { value };
}

inline void
Context::set(const char *parameter, std::vector<double> values)
{
  parameters[parameter] = values;
}

inline void
Context::set(const std::string& parameter, std::vector<double> values)
{
  parameters[parameter] = values;
}

inline void
Context::set_if_not_set(const char *key, std::vector<double> values)
{
  if (parameters.find(key) == parameters.end())
    set(key, values);
}

inline void
Context::push_back(const std::string& parameter, double value)
{
  parameters[parameter].push_back(value);
}

inline size_t
Context::size(const char *parameter)
{
  return parameters.at(parameter).size();
}

inline void
Context::set_time_adjust(const std::string & parameter,
		double time_period,
		size_t from,
		size_t step,
		TimeAdjustMethod method)
{
  parameters_to_time_adjust[parameter] =
    std::tuple<double,size_t, size_t, TimeAdjustMethod>
    (time_period, from, step, method);
}

inline double
Context::get_time_adjust(std::string & s)
{
  return std::get<0>(parameters_to_time_adjust[s]);
}


// Reporter inline methods


inline void
Reporter::setSize(const size_t numVals)
{
  values_.resize(numVals);
}

inline void
Reporter::getReportValue(const Simulation &s, size_t id)
{
  values_[id] = getValueFunc_(s);
}

inline void
Reporter:: calculate()
{
  calcVal_ = calcFunc_(values_);
}

inline void
Reporter::output()
{
  outputFunc_(calcVal_);
}


///////////////////

// Options inline methods

inline Options::Options()
  : events_({advanceTimeEvent})
  , tests_({})
  , argc_(0)
  , argv_(NULL)
  , individual_agent_create_func_(create_default_agent)
  , all_agents_create_func_(NULL)
  , before_all_simulations_func_(NULL)
  , before_each_simulation_func_(NULL)
  , after_each_simulation_func_(NULL)
{}
inline Options&
Options::events(Events events)
{ events_ = events; return *this; }

inline Options&
Options::additionalTests(Tests tests)
{ tests_ = tests; return *this; }

inline Options&
Options::commandLine(int argc, char **argv)
{ argc_ = argc; argv_ = argv; return *this; }

inline Options&
Options::agentCreate(std::function<Agent *(Context &)>
		     individual_agent_create_func)
{ individual_agent_create_func_ = individual_agent_create_func; return *this; }
inline Options&
Options::allAgentsCreate(std::function<void(Simulation &)>
			 all_agents_create_func)
{ all_agents_create_func_ = all_agents_create_func; return *this; }
inline Options&
Options::beforeAllSimulations(std::function<void(Simulation &)>
			      before_all_simulations_func)
{
  before_all_simulations_func_ = before_all_simulations_func;
  return *this;
}
inline Options&
Options::beforeEachSimulation(std::function<void(Simulation &)>
			     before_each_simulation_func)
{
  before_each_simulation_func_ = before_each_simulation_func;
  return *this;
}
inline Options&
Options::afterEachSimulation(std::function<void(Simulation &)>
			     after_each_simulation_func)
{
  after_each_simulation_func_ = after_each_simulation_func;
  return *this;
}
inline Options&
Options::parameter(const char *parameter,
		   const std::vector<double> & values)
{
  context_.set(parameter, values);
  return *this;
}
inline Options&
Options::timeAdjust(const std::string & parameter,
		    double time_period,
		    size_t from,
		    size_t step,
		    TimeAdjustMethod method)
{
  context_.set_time_adjust(parameter, time_period, from, step, method);
  return *this;
}
inline Options&
Options::report(std::function<double(const Simulation &)> getValueFunc,
		std::function<double(std::vector<double> &)> calcFunc,
		std::function<void(const double)> outputFunc)
{
  reporters_.push_back(Reporter(getValueFunc, calcFunc, outputFunc));
  return *this;
}

///////////////////

// Simulation inline methods

inline
Simulation::Simulation() {}

inline
Simulation::Simulation(const Options & options)
{
  setOptions(options);
}

inline void
Simulation::init(unsigned seed,
		 Events events,
		 std::function<Agent *(Context &)> agent_create_func,
		 std::function<void (Simulation &)> create_all_agents_func)
{
  events_ = events;
  agent_create_func_ = agent_create_func;
  create_all_agents_func_ = create_all_agents_func;
  init(seed);
}

inline double
Simulation::time_correct_prob(const double parameter_prob,
			      const double parameter_time_period)
{
  return sim::time_correct_prob(parameter_prob,
				parameter_time_period,
				time_step);
}

inline void
Simulation::move_agent(std::vector <Agent *> &from_agents,
		       std::vector <Agent *> &to_agents,
		       const std::size_t agent_index)
{
  to_agents.push_back(from_agents[agent_index]);
  from_agents[agent_index] = from_agents.back();
  from_agents.pop_back();
}

inline void
Simulation::remove_dead_agents()
{
  for (size_t i = 0; i < agents.size(); ++i)
    if (agents[i]->alive == false)
      move_agent(agents, dead_agents, i);
}

inline
void Simulation::shuffle_agents()
{
  shuffle(agents.begin(), agents.end(), rng);
}


#endif
