# include "sim.hh"

using namespace sim;

namespace sim {
  thread_local  std::mt19937_64 rng;
  std::vector< Tracking > trackings;
}

//////////////////////////////
/* COMMAND LINE PROCESSING */

/*
  Simple command line processing functions taken from:
  http://stackoverflow.com/questions/865668/parse-command-line-arguments
*/

char* getCmdOption(char ** begin, char ** end, const std::string & option)
{
  char ** itr = std::find(begin, end, option);
  if (itr != end && ++itr != end)
    {
      return *itr;
    }
  return 0;
}

bool cmdOptionExists(char** begin, char** end, const std::string& option)
{
  return std::find(begin, end, option) != end;
}

/* END COMMAND LINE PROCESSING */
/////////////////////////////////


bool
sim::is_event(const double risk)
{
  std::uniform_real_distribution<double> uni;
  if (uni(rng) < risk)
    return true;
  else
    return false;
}

void
sim::advanceTimeEvent(Simulation &simulation)
{
  simulation.current_date +=  simulation.time_step;
}

void
sim::deathEvent(Simulation &simulation) {}


std::vector< std::vector<std::string> >
get_csv_table(const char *parameter_file_name)
{
  std::vector< std::vector<std::string> > csv_table;
  std::ifstream infile;

  std::string line, cell;

  infile.open(parameter_file_name);
  if ( infile.is_open() ) {
    while ( std::getline(infile, line) ) {
      std::istringstream line_stream(line);
      std::vector<std::string> cells;
      while ( std::getline(line_stream, cell, ',') )
	if (cell[0] == '#')
	  break;
	else if (cell.size())
	  cells.push_back(cell);
      csv_table.push_back(cells);
    }
    infile.close();
  } else {
    std::cout << "Could not open file "
	      << parameter_file_name << "." << std::endl;
  }
  return csv_table;
}

void
process_parameter_file(Context& context,
		       const char * parameter_file_name)
{
  std::vector< std::vector<std::string> > csv_table =
    get_csv_table(parameter_file_name);

  for (auto line : csv_table) {
    std::vector<double> values;
    if ( line.size() > 1 ) {
      if (line[0] == "@") {
	size_t num_agents = std::stoul(line[1]);
	std::vector<std::string> characteristics;
	for (auto it = line.begin() + 2; it != line.end(); ++it) {
	  characteristics.push_back(*it);
	}
	context.set_initial(num_agents, characteristics);
      } else {
	for (auto it = line.begin() + 1; it != line.end(); ++it) {
	  if ((*it)[0] == '!' || (*it)[0] == '^' || (*it)[0] == '&') {
	    // 01234567
	    // !1.0!0!2
	    double time_period = 1.0;
	    size_t from = 0, step = 1;
	    TimeAdjustMethod method;
	    switch ((*it)[0]) {
	    case '!': method = LINEAR; break;
	    case '&' : method = PROBABILITY; break;
	    case '^' : method = COMPOUND; break;
	    }
	    size_t index1 = (*it).find('!', 1), index2 = std::string::npos;
	    if (index1 != std::string::npos) {
	      time_period = std::stod((*it).substr(1, index1 - 1));
	      index2 = (*it).find('!', index1 + 1);
	      if (index2 != std::string::npos){
		from = std::stoul((*it).substr(index1 + 1, index2 - index1));
		step = std::stoul((*it).substr(index2 + 1));
	      } else {
		from = std::stod((*it).substr(index1 + 1));
	      }
	    } else {
	      time_period = std::stod((*it).substr(1));
	    }
	    context.set_time_adjust(line[0], time_period, from, step, method);
	  } else {
	    values.push_back(std::stod(*it));
	  }
	}
	context.set(line[0], values);
      }
    }
  }
}

void
process_command_line_parameters(Context& context, const char *parameter_str)
{
  std::vector<std::string> parm_line;
  std::string s(parameter_str);
  size_t i = 0, j = 0;

  do {
    j = s.find(":", i);
    if (j == std::string::npos)
      j = s.size();
    std::string p = s.substr(i, j - i);
    i = j + 1;
    parm_line.push_back(p);
  } while(j != s.size());

  for (auto p : parm_line) {
    if ( (p[0] >= 'a' && p[0] <= 'z') ||
	 (p[0] >= 'A' && p[0] <= 'Z')) {
      s = p;
      context.set(s, {});
    } else if (p[0] == '!') {
      context.set_time_adjust(s, stod(p.substr(1)));
    } else {
      context.push_back(s, stod(p));
    }
  }
}

class Report {
public:
  Report(tst::TestSeries &test_series,
	 size_t num_agents) : t_(test_series), n_(num_agents) {};
  void operator()(sim::Simulation &s)
  {
    TESTEQ(t_, s.agents.size(), n_, "number of agents after simulation");
    TESTDBL(t_, s.current_date,
	    (s.context("START_DATE") +
	     (s.context("ITERATIONS") * s.context("TIME_STEP"))),
	    "end and start dates correspond to time steps");
  }
private:
  tst::TestSeries &t_;
  size_t n_;
};

void run_tests(const Tests & user_tests)
{
  Simulation s;
  tst::TestSeries t;

  // Test auxiliary functions

  {
    std::vector<double> x = {};
    TESTDBL(t, median(x), 0.0, "median of 0");
    std::vector<double> a = {4.0, 3.0, 6.0};
    TESTDBL(t, median(a), 4.0, "median of 3");
    std::vector<double> b = {4.0, 3.0, 6.0, 5.0};
    TESTDBL(t, median(b), 4.5, "median of 4");
    std::vector<double> c = {3,4,4,5,5,4,3,7,6,9,6,8};
    TESTDBL(t, median(c), 5.0, "median of 12");
  }

  // Test simulation
  process_parameter_file(s.context, "test_csv.csv");
  s.context.set_defaults_not_yet_set();

  TESTDBL(t, s.context("NEW_PARM", 0), 5.0, "parameter type added");
  TESTDBL(t, s.context("NEW_PARM", 1), 1.0, "multiple entries added");
  TESTDBL(t, s.context("HIV_PREVALENCE_STAGE",0), 0.2, "entry correct 1");
  TESTDBL(t, s.context("HIV_PREVALENCE_STAGE",1), 0.1, "entry correct 2");
  TESTDBL(t, s.context("HIV_PREVALENCE_STAGE",2), 0.3, "entry correct 3");
  TESTDBL(t, s.context("HIV_PREVALENCE_STAGE",3), 0.2, "entry correct 4");
  TESTDBL(t, s.context.get_time_adjust("HIV_PREVALENCE_STAGE").time_period, 1.0,
	  "time adjuster set");
  TESTDBL(t, s.context.get_time_adjust("HIV_PREVALENCE_STAGE").from, 0,
	  "from set");
  TESTDBL(t, s.context.get_time_adjust("HIV_PREVALENCE_STAGE").step, 1,
	  "step set");
  TESTDBL(t, s.context.get_time_adjust("HIV_RISK_DEATH").time_period, 1.25,
	  "time adjuster set");
  TESTDBL(t, s.context.get_time_adjust("HIV_RISK_DEATH").from, 1,
	  "from set");
  TESTDBL(t, s.context.get_time_adjust("HIV_RISK_DEATH").step, 2,
	  "step set");
  s.context.adjust_parameters_to_time_period();
  TESTDBL(t, s.context("HIV_PREVALENCE_STAGE", 1),
	  0.00202411247526512738659221213311, "time corrected probability");
  TESTDBL(t, s.context("HIV_RISK_DEATH", 0), 500,
	  "index into time adjuster left alone 1");
  TESTDBL(t, s.context("HIV_RISK_DEATH", 2), 400,
	  "index into time adjuster left alone 2");
  TESTDBL(t, s.context("HIV_RISK_DEATH", 4), 350,
	  "index into time adjuster left alone 3");
  TESTDBL(t, s.context("HIV_RISK_DEATH", 6), 200,
	  "index into time adjuster left alone 4");
  TESTDBL(t, s.context("HIV_RISK_DEATH", 8), 50,
	  "index into time adjuster left alone 4");
  TESTDBL(t, s.context("HIV_RISK_DEATH", 9), 0.0139978861008127619669494379195,
	  "correct time adjustment");

  try {
    s.context("HIV_PREVALENCE_STAGE",4);
    TESTEQ(t, 0, 1, "exception not thrown 1");
  } catch (InvalidParameter &e) {
    TESTEQ(t, 1, 1, "exception thrown 1");
  }
  try {
    s.context("HIV_PREVALENCE_");
    TESTEQ(t, 0, 1, "exception not thrown 2");
  } catch(InvalidParameter &e) {
    TESTEQ(t, 1, 1, "exception thrown 2");
  }

  size_t n = s.context("NUM_AGENTS");
  s.setOptions(Options()
	       .events({advanceTimeEvent, deathEvent})
	       .afterEachSimulation(Report(t, n)));
  s.simulate();

  for (auto test_func : user_tests)
    test_func(t);

  t.summary();
}


void
sim::process_command_line(Context &context,
			  int argc,
			  char *argv[],
			  const Tests& user_tests)
{
  // Command line options
  bool test = cmdOptionExists(argv, argv + argc, "-t");
  bool test_only = cmdOptionExists(argv, argv + argc, "-to");
  bool verbose = cmdOptionExists(argv, argv + argc, "-v");
  char *parameter_file_str = getCmdOption(argv, argv + argc, "-f");
  char *verbosity_str = getCmdOption(argv, argv + argc, "-v");
  char *parameter_str = getCmdOption(argv, argv + argc, "-p");

  if (test || test_only)
    run_tests(user_tests);

  if (test_only)
    context.set("_SIM", 0.0);
  else
    context.set("_SIM", 1.0);

  if (parameter_file_str)
    process_parameter_file(context, parameter_file_str);

  if (parameter_str)
    process_command_line_parameters(context, parameter_str);

  // Determine whether or not to print the parameters
  bool print_parms = false;
  if (verbose) {
    if (verbosity_str) {
      std::string s(verbosity_str);
      if (s.find('p') != std::string::npos) {
	print_parms = true;
      };
    } else {
      print_parms = true;
    }
    if (print_parms) {
      std::cout << "************************" << std::endl;
      std::cout << "Parameters" << std::endl;
      context.print_parameters();
      std::cout << "************************" << std::endl;
    }
  }
}

Agent::Agent(Context &c) : context_(c)
{
  std::uniform_real_distribution<double> uni;
  id = context_.last_agent++;
  alive = true;
}

double
sim::time_correct_linear(const double parameter_prob,
			 const double parameter_time_period,
			 const double actual_time_period)
{
  return parameter_prob * actual_time_period / parameter_time_period;
}

double
sim::time_correct_compound(const double parameter_value,
			   const double parameter_time_period,
			   const double actual_time_period)
{
  double k = actual_time_period / parameter_time_period;
  double t = pow(parameter_value, k);
  return t;
}


double
sim::time_correct_prob(const double parameter_prob,
		       const double parameter_time_period,
		       const double actual_time_period)
{

  double t = actual_time_period / parameter_time_period;
  double not_p = 1 - parameter_prob;
  double not_p_corrected = pow(not_p, t);
  double new_p = 1 - not_p_corrected;
  return new_p;
}


double sim::Agent::age(Simulation &s)
{
  return s.current_date - dob;
}

void sim::Agent::die(Simulation &s, const std::string & c)
{
  alive = false;
  dod = s.current_date;
  cause = c;
}

Agent* sim::create_default_agent(Context & c)
{
  std::uniform_real_distribution<double> uni;
  Agent *a  = new Agent(c);

  a->sex = uni(rng) < c("PROB_MALE") ?
		      MALE : FEMALE;
  { // dob
    std::uniform_real_distribution<double>
      uni_age(c("EARLIEST_BIRTH_DATE"),
	      c("LATEST_BIRTH_DATE"));
    a->dob = uni_age(rng);
  }
  return a;
}


double sim::numAgents(const Simulation &s)
{
  return (double) s.agents.size();
}

// Context methods

void
Context::set_defaults_not_yet_set()
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
  set_if_not_set("SIMULATIONS_PER_THREAD", { 10.0 });
  set_if_not_set("THREADED", { 1.0 });
}

double
Context::operator()(const char * key, size_t i) const
{
  if (parameters.find(key) == parameters.end()) {
    std::stringstream ss;
#ifdef __GNUC__
#ifdef TRACING
    void *array[20];
    int size;
    size = backtrace(array, 20);
    backtrace_symbols_fd(array, size, STDERR_FILENO);
#endif
#endif
    ss << "Unknown parameter: " << key << std::endl;
    throw InvalidParameter(ss.str());
  } else {
    auto p = parameters.at(key);
    if (i >= p.size())  {
      std::stringstream ss;
      ss << "Index " << i << " out of range for parameter "
	 << key << std::endl;
      throw  InvalidParameter(ss.str());
    } else {
      return p[i];
    }
  }
}

std::vector<double> &
Context::get(const char *key)
{
  auto v = parameters.find(key);
  if (v == parameters.end()) {
    std::stringstream ss;
    ss << "Unknown parameter: " << key << std::endl;
    throw InvalidParameter(ss.str());
  }
  return v->second;
}


TimeAdjust
Context::get_time_adjust(const char * s)
{
  TimeAdjust output;
  auto t = parameters_to_time_adjust[s];
  output.time_period = std::get<0>(t);
  output.from = std::get<1>(t);
  output.step = std::get<2>(t);
  output.method = std::get<3>(t);
  return output;
}


void
Context::set_initial(size_t num_agents,
		     std::vector< std::string > characteristics)
{
  initialization.push_back(std::pair<size_t, std::vector<std::string> >
			   (num_agents, characteristics));
}

std::pair<size_t, std::vector<std::string> > &
Context::get_initial(size_t i)
{
  return initialization[i];
}

std::vector< std::pair<size_t, std::vector<std::string> > > &
Context::get_initial()
{
  return initialization;
}

void
Context::convert_probabilities_to_time_period(const char * key,
					      double parameter_time_period,
					      size_t start,
					      size_t step)
{
  double actual_time_period = (*this)("TIME_STEP");
  std::vector<double> & values = get(key);
  for (auto it = values.begin() + start; it < values.end(); it += step)
    *it = time_correct_prob(*it, parameter_time_period, actual_time_period);
}


void
Context::convert_linear_values_to_time_period(const char * key,
					     double parameter_time_period,
					     size_t start,
					     size_t step)
{
  double actual_time_period = (*this)("TIME_STEP");
  std::vector<double> & values = get(key);
  for (auto it = values.begin() + start; it < values.end(); it += step)
    *it = time_correct_linear(*it, parameter_time_period,
			      actual_time_period);
}

void
Context::convert_compound_product_to_time_period(const char * key,
						 double parameter_time_period,
						 size_t start,
						 size_t step)
{
  double actual_time_period = (*this)("TIME_STEP");
  std::vector<double> & values = get(key);
  for (auto it = values.begin() + start; it < values.end(); it += step)
    *it = time_correct_compound(*it, parameter_time_period,
				actual_time_period);
}

void
Context::adjust_parameters_to_time_period()
{
  for (auto entry : parameters_to_time_adjust) {
    if (std::get<3>(entry.second) == PROBABILITY)
      convert_probabilities_to_time_period(entry.first.c_str(),
					   std::get<0>(entry.second),
					   std::get<1>(entry.second),
					   std::get<2>(entry.second));
    else if (std::get<3>(entry.second) == LINEAR)
      convert_linear_values_to_time_period(entry.first.c_str(),
					   std::get<0>(entry.second),
					   std::get<1>(entry.second),
					   std::get<2>(entry.second));
    else
      convert_compound_product_to_time_period(entry.first.c_str(),
					      std::get<0>(entry.second),
					      std::get<1>(entry.second),
					      std::get<2>(entry.second));
  }
}

void
Context::print_parameters()
{
  for (auto entry : parameters) {
    std::cout << entry.first << ":\t";
    for (auto d : entry.second)
      std::cout << d << "\t";
    std::cout << std::endl;
  }
}

///////////////


// Reporter methods

Reporter::Reporter(std::function<double(const Simulation &)> getValueFunc,
		   std::function<double(std::vector<double> &)> calcFunc,
		   std::function<void(const double)> outputFunc)
{
  getValueFunc_ = getValueFunc;
  calcFunc_ = calcFunc;
  outputFunc_ = outputFunc;
}


//////////////


// Simulation methods

Simulation::~Simulation()
{
  for (auto a : agents)
    delete a;
  for (auto a : dead_agents)
    delete a;
  for (auto a : aged_out_agents)
    delete a;
  if (simulation_num == 0 && reporters_) {
    delete reporters_;
    reporters_ = NULL;
  }
}

void
Simulation::setOptions(const Options & options)
{
  events_ = options.events_;
  tests_ = options.tests_;
  argc_ = options.argc_;
  argv_ = options.argv_;
  tracking_on_ = options.tracking_on_;
  agent_create_func_ = options.individual_agent_create_func_;
  create_all_agents_func_ = options.all_agents_create_func_;
  before_all_simulations_func_ = options.before_all_simulations_func_;
  before_each_simulation_func_ = options.before_each_simulation_func_;
  after_each_simulation_func_ = options.after_each_simulation_func_;
  context = options.context_;
  if (options.reporters_.size()) {
    reporters_ = new std::vector<Reporter>;
    for (auto & r : options.reporters_)
      reporters_->push_back(r);
  }
}


void
Simulation::init(unsigned seed)
{
  rng.seed(seed);
  time_step = context("TIME_STEP");
  current_date = context("START_DATE");
  total_iterations = context("ITERATIONS");
  if (create_all_agents_func_)
    create_all_agents_func_(*this);
  else
    for (unsigned i = 0; i < context("NUM_AGENTS"); ++i) {
      Agent *a = agent_create_func_(context);
      agents.push_back(a);
    }
}







void
Simulation::simulate()
{
  std::vector<std::thread> threads;
  size_t num_sims, sim_per_thread, num_threads;

  if (argc_) {
    process_command_line(context, argc_, argv_, tests_);
    if (context("_SIM") == 0.0)
      return;
  }

  context.set_defaults_not_yet_set();
  context.adjust_parameters_to_time_period();

  num_sims = context("NUM_SIMULATIONS");

  if (tracking_on_)
    trackings.reserve(num_sims);

  if (before_all_simulations_func_)
    before_all_simulations_func_(*this);

  if (reporters_)
    for (auto & r : *reporters_)
      r.setSize(num_sims);

  if (context("THREADED")) {
    sim_per_thread = context("SIMULATIONS_PER_THREAD");
    num_threads = ceil((double) num_sims / sim_per_thread );
    for (size_t i = 0; i < num_threads; ++i) {
      size_t from = i * sim_per_thread;
      size_t to = std::min((i + 1) * sim_per_thread, num_sims);
      threads.push_back(std::thread(threaded_part_of_sim_loop,
				    std::ref(*this), from, to));
    }
    for (auto & t : threads)
      t.join();
  } else {
    threaded_part_of_sim_loop(*this, 0, num_sims);
  }
  if (reporters_)
    for (auto & r : *reporters_) {
      r.calculate();
      r.output();
    }
}


void
Simulation::simulate_once()
{
  for (current_iteration = 0; current_iteration < total_iterations;
       ++current_iteration)
    for (auto event : events_)
      event(*this);
}
