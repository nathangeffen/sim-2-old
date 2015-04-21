#ifndef __SIM_H__
#define __SIM_H__

#include <cassert>
#include <cfloat>
#include <cmath>

#include <array>
#include <algorithm>
#include <exception>
#include <iomanip>
#include <iostream>
#include <fstream>
#include <functional>
#include <random>
#include <sstream>
#include <string>
#include <vector>
#include <unordered_map>

namespace sim {

  class Simulation;

  const double YEAR = 1.0;
  const double MONTH = 1.0 / 12.0;
  const double WEEK = 1.0 / 52.0;
  const double DAY = 1.0 / 365.25;

  extern std::mt19937 rng;

  bool is_event(const double risk);

  typedef std::vector< std::function<void(Simulation &)> > Events;

  double
  time_correct_prob(const double parameter_prob,
		    const double parameter_time_period,
		    const double actual_time_period);
  double
  time_correct_prob(const double parameter_prob,
		    const double parameter_time_period);


  enum Sex {
    MALE = 0,
    FEMALE = 1
  };

  struct InvalidParameter : public std::exception
  {
  public:
    InvalidParameter(const std::string & s =
		     "Invalid index into parameter array")
    {
      msg = s;
    }
  private:
    std::string msg;
    const char * what () const throw ()
    {
      return msg.c_str();
    }
  };

  struct TimeAdjust {
    double time_period;
    size_t from;
    size_t step;
  };

  class Common {
  public:
    void set_defaults_not_yet_set()
    {
      set_if_not_set("NUM_AGENTS", { 100 });
      set_if_not_set("TIME_STEP", { DAY });
      set_if_not_set("NUM_YEARS", { 20.0 });
      set_if_not_set("ITERATIONS", { (*this)("NUM_YEARS") *
	    (1.0 /  (*this)("TIME_STEP")) });
      set_if_not_set("PROB_MALE", { 0.5 });
      set_if_not_set("START_DATE", { 2010.0 });
      set_if_not_set("YOUNGEST_AGE", { 15.0 });
      set_if_not_set("OLDEST_AGE", { 49.0 });
      set_if_not_set("LATEST_BIRTH_DATE",
		     { (*this)("START_DATE") - (*this)("YOUNGEST_AGE") });
      set_if_not_set("EARLIEST_BIRTH_DATE",
		     { (*this)("START_DATE") - (*this)("OLDEST_AGE") });
      set_if_not_set("HIV_PREVALENCE_STAGE", { 0.025, 0.025, 0.025, 0.025 });
      set_if_not_set("HIV_PREVALENCE", { 0.19 });
      set_if_not_set("PROB_CIRCUMCISED", { 0.2 });
      set_if_not_set("NUM_SIMULATIONS", { 1.0 });
    }
    double operator()(const char * c, size_t i)
    {
      try {
	auto p = parameters.at(c);
	if (i >= p.size())  {
	  std::stringstream ss;
	  ss << "Index " << i << " out of range for parameter "
	     << c << std::endl;
	  throw  InvalidParameter(ss.str());
	} else {
	  return p[i];
	}
      } catch (const std::out_of_range& e) {
	std::stringstream ss;
	ss << "Unknown parameter: " << c << std::endl;
	throw InvalidParameter(ss.str());
      }
    }
    double operator()(const char * c)
    {
      return operator()(c, 0);
    }
    double operator()(const std::string & s) {
      return operator()(s, 0);
    }
    double operator()(const std::string &s, size_t i)
    {
      return operator()(s.c_str(), i);
    }
    void set(const char *parameter, double value)
    {
      parameters[parameter] = { value };
    }
    void set(const char *parameter, std::vector<double> values)
    {
      parameters[parameter] = values;
    }
    void set(const std::string& parameter, std::vector<double> values)
    {
      parameters[parameter] = values;
    }
    void
    set_if_not_set(const char *key, std::vector<double> values)
    {
      if (parameters.find(key) == parameters.end())
	set(key, values);
    }
    std::vector<double> & get(const char *key)
    {
      auto v = parameters.find(key);
      if (v == parameters.end()) {
	std::stringstream ss;
	ss << "Unknown parameter: " << key << std::endl;
	throw InvalidParameter(ss.str());
      }
      return v->second;
    }
    void push_back(const std::string& parameter, double value)
    {
      parameters[parameter].push_back(value);
    }
    size_t size(const char *parameter)
    {
      return parameters.at(parameter).size();
    }
    void
    set_time_adjust(const std::string & parameter,
		    double time_period,
		    size_t from = 0.0,
		    size_t step = 0.0)
    {
      parameters_to_time_adjust[parameter] =
	std::tuple<double,size_t, size_t>(time_period, from, step);
    }
    TimeAdjust get_time_adjust(const char * s)
    {
      TimeAdjust output;
      auto t = parameters_to_time_adjust[s];
      output.time_period = std::get<0>(t);
      output.from = std::get<1>(t);
      output.step = std::get<2>(t);
      return output;
    }
    double get_time_adjust(std::string & s)
    {
      return std::get<0>(parameters_to_time_adjust[s]);
    }
    void
    convert_probabilities_to_time_period(const char * key,
					 double parameter_time_period,
					 size_t start = 0,
					 size_t step = 1)
    {
      double actual_time_period = (*this)("TIME_STEP");
      std::vector<double> & values = get(key);
      for (auto it = values.begin() + start; it < values.end(); it += step)
	*it = time_correct_prob(*it, parameter_time_period, actual_time_period);
    }

    void adjust_parameters_to_time_period()
    {
      for (auto entry : parameters_to_time_adjust)
	convert_probabilities_to_time_period(entry.first.c_str(),
					     std::get<0>(entry.second),
					     std::get<1>(entry.second),
					     std::get<2>(entry.second));
    }

    void print_parameters()
    {
      for (auto entry : parameters) {
	std::cout << entry.first << ":\t";
	for (auto d : entry.second)
	  std::cout << d << "\t";
	std::cout << std::endl;
      }
    }
    unsigned last_agent = 0;
  private:
    std::unordered_map<std::string, std::vector<double> > parameters;
    std::unordered_map<std::string, std::tuple<double, size_t, size_t> >
    parameters_to_time_adjust;
  };

  class Agent
  {
  public:
    Agent(Common &common,
	  std::function<void(Agent *, Common &)> agent_initiation_func)
      : common_(common) {
      agent_initiation_func(this, common_);
    }
    double age(Simulation &s);
    void die(Simulation &s, const std::string& cause);
    unsigned id;
    bool alive;
    double dob;
    double dod;
    std::string cause;
    Sex sex;
    int hiv;
    double cd4;
    bool circumcised;
    double orientation;
    double riskiness;
    unsigned num_partners;
    std::vector<unsigned> partners;
  private:
    Common &common_;
  };

  void default_agent_initiation(Agent *a, Common &c);

  class Simulation
  {
  public:
    Simulation(Common &c) : common(c) {}
    ~Simulation()
    {
      for (auto a : agents)
	delete a;
      for (auto a : dead_agents)
	delete a;
      for (auto a : aged_out_agents)
	delete a;
    }
    void init(Events e,
	      std::function<void(Agent *, Common &)> agent_initiation_func =
	      default_agent_initiation)
    {
      time_step = common("TIME_STEP");
      current_date = common("START_DATE");
      events = e;
      for (unsigned i = 0; i < common("NUM_AGENTS"); ++i) {
	Agent *a = new Agent(common, agent_initiation_func);
	agents.push_back(a);
      }
    }

    double
    time_correct_prob(const double parameter_prob,
		      const double parameter_time_period)
    {
      return sim::time_correct_prob(parameter_prob,
				    parameter_time_period,
				    time_step);
    }

    void
    move_agent(std::vector <Agent *> &from_agents,
	       std::vector <Agent *> &to_agents,
	       const std::size_t agent_index)
    {
      to_agents.push_back(from_agents[agent_index]);
      from_agents[agent_index] = from_agents.back();
      from_agents.pop_back();
    }

    void remove_dead_agents()
    {
      for (size_t i = 0; i < agents.size(); ++i)
	if (agents[i]->alive == false)
	  move_agent(agents, dead_agents, i);
    }

    void simulate()
    {
      unsigned iterations = common("ITERATIONS");
      for (iteration = 0; iteration < iterations; ++iteration)
	for (auto event : events)
	  event(*this);
    }
    std::vector<Agent *> agents;
    std::vector<Agent *> aged_out_agents;
    std::vector<Agent *> dead_agents;
    Events events;
    double current_date, time_step;
    unsigned iteration;
    Common &common;
  };
  void process_command_line(Common &common, int argc, char *argv[]);
  void advanceTimeEvent(Simulation &simulation);
  void deathEvent(Simulation &simulation);
}

#endif
