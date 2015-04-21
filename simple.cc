#include "sim.hh"


void set_parameters(sim::Common &c)
{
  c.set("AVG_CD4_DECLINE", 67.2);
  c.set("RISK_DEATH_MALE",
	{
	    0.000436,
	    0.000304,
	    0.000232,
	    0.000172,
	    0.000155,
	    0.000143,
	    0.000131,
	    0.000115,
	    0.000096,
	    0.000082,
	    0.000086,
	    0.000125,
	    0.000205,
	    0.000319,
	    0.000441,
	    0.000562,
	    0.00069,
	    0.00082,
	    0.000949,
	    0.001085,
	    0.001213,
	    0.001304,
	    0.001345,
	    0.00135,
	    0.001342,
	    0.00134,
	    0.001342,
	    0.001356,
	    0.00138,
	    0.001408,
	    0.001435,
	    0.001466,
	    0.001499,
	    0.001539,
	    0.001592,
	    0.00166,
	    0.001741,
	    0.001837,
	    0.001953,
	    0.002084,
	    0.002241,
	    0.002439,
	    0.002686,
	    0.002975,
	    0.003297,
	    0.003639,
	    0.003997,
	    0.004366,
	    0.00475,
	    0.005156,
	    0.005596,
	    0.006078,
	    0.006605,
	    0.007174,
	    0.007805,
	    0.008464,
	    0.009095,
	    0.009676,
	    0.010245,
	    0.010865,
	    0.011592,
	    0.012444,
	    0.013451,
	    0.014608,
	    0.015927,
	    0.01737,
	    0.018895,
	    0.020484,
	    0.022191,
	    0.024139,
	    0.026364,
	    0.028808,
	    0.03148,
	    0.034442,
	    0.037855,
	    0.041725,
	    0.045932,
	    0.050469,
	    0.055465,
	    0.061179,
	    0.067698,
	    0.074923,
	    0.082891,
	    0.091725,
	    0.101575,
	    0.112568,
	    0.124795,
	    0.138305,
	    0.153107,
	    0.169195,
	    0.186543,
	    0.205115,
	    0.224867,
	    0.245744,
	    0.266454,
	    0.286625,
	    0.305869,
	    0.323783,
	    0.339972,
	    0.356971,
	    0.374819,
	    0.39356,
	    0.413238,
	    0.4339,
	    0.455595,
	    0.478375,
	    0.502293,
	    0.527408,
	    0.553778,
	    0.581467,
	    0.610541,
	    0.641068,
	    0.673121,
	    0.706777,
	    0.742116,
	    0.779222,
	    0.818183,
	    0.859092,
	    0.902047
	    }
	);
  c.set("RISK_DEATH_FEMALE",
	{
	  0.005562,
	    0.000396,
	    0.000214,
	    0.000162,
	    0.000132,
	    0.000117,
	    0.000106,
	    0.000099,
	    0.000093,
	    0.00009,
	    0.00009,
	    0.000096,
	    0.000111,
	    0.000137,
	    0.00017,
	    0.000207,
	    0.000245,
	    0.000282,
	    0.000318,
	    0.000352,
	    0.000388,
	    0.000423,
	    0.000454,
	    0.000476,
	    0.000494,
	    0.000511,
	    0.000531,
	    0.000553,
	    0.000579,
	    0.000608,
	    0.000641,
	    0.000677,
	    0.000719,
	    0.000765,
	    0.000818,
	    0.000879,
	    0.000948,
	    0.001022,
	    0.0011,
	    0.001185,
	    0.001279,
	    0.001387,
	    0.001518,
	    0.001676,
	    0.001858,
	    0.002055,
	    0.002262,
	    0.00248,
	    0.002709,
	    0.002947,
	    0.003209,
	    0.003484,
	    0.003751,
	    0.004,
	    0.004246,
	    0.00452,
	    0.004836,
	    0.005185,
	    0.00557,
	    0.006001,
	    0.006489,
	    0.007046,
	    0.007686,
	    0.008419,
	    0.009249,
	    0.010201,
	    0.011255,
	    0.012372,
	    0.013538,
	    0.014793,
	    0.016233,
	    0.017882,
	    0.019693,
	    0.021671,
	    0.023866,
	    0.026437,
	    0.029368,
	    0.032519,
	    0.03587,
	    0.039555,
	    0.043828,
	    0.048808,
	    0.054434,
	    0.060762,
	    0.067889,
	    0.075926,
	    0.084968,
	    0.095093,
	    0.106352,
	    0.118777,
	    0.132384,
	    0.147181,
	    0.163161,
	    0.180314,
	    0.198615,
	    0.217125,
	    0.235558,
	    0.253602,
	    0.270923,
	    0.287178,
	    0.304409,
	    0.322673,
	    0.342033,
	    0.362555,
	    0.384309,
	    0.407367,
	    0.431809,
	    0.457718,
	    0.485181,
	    0.514292,
	    0.545149,
	    0.577858,
	    0.61253,
	    0.649282,
	    0.688238,
	    0.729533,
	    0.773305,
	    0.818183,
	    0.859092,
	    0.902047
	});

  c.set("HIV_RISK_DEATH",
	{
	  500.0,
	    0.01,
	    400.0,
	    0.02,
	    350.0,
	    0.03,
	    200.0,
	    0.4,
	    50.0,
	    0.6
	});
}


void cd4Event(sim::Simulation &s)
{
  static thread_local double avg_cd4_decline = s.common("AVG_CD4_DECLINE");
  for (auto agent : s.agents) {
    if (agent->hiv) {
      agent->cd4 = std::max(0.0,
			    agent->cd4 -
			    avg_cd4_decline * s.time_step);
    }
  }
}


void deathEvent(sim::Simulation &s)
{
  static thread_local std::vector<double> & male_mort =
    s.common.get("RISK_DEATH_MALE");
  static thread_local std::vector<double> & female_mort=
    s.common.get("RISK_DEATH_FEMALE");
  static thread_local std::vector<double> & hiv_mort =
    s.common.get("HIV_RISK_DEATH");
  static thread_local size_t size = hiv_mort.size();

  for (auto agent : s.agents) {
    unsigned age = agent->age(s);
    double risk_of_death;

    // Risk of death for everyone due to age
    if (agent->sex == sim::MALE)
      risk_of_death = male_mort[age];
    else
      risk_of_death = female_mort[age];
    if (sim::is_event(risk_of_death))
      agent->die(s, "AGE");

    // HIV risk of death
    if (agent->alive && agent->hiv) {
      for (size_t i = size; i > 0; i -= 2) {
	if (agent->cd4 < hiv_mort[i - 2]) {
	  double hiv_risk_of_death = hiv_mort[i - 1];
	  if (sim::is_event(hiv_risk_of_death)) {
	    agent->die(s, "HIV");
	    break;
	  }
	}
      }
    }
  }
  s.remove_dead_agents();
}


void report(sim::Simulation &s, unsigned sim_number)
{
  double prop_dead = (double) s.dead_agents.size() /
    (s.agents.size() + s.dead_agents.size());
  std::cout << sim_number << ", alive, " << s.agents.size() << std::endl;
  std::cout << sim_number << ", dead, " << s.dead_agents.size() << std::endl;
  std::cout << sim_number << ", proportion dead, " << prop_dead << std::endl;
}

int main(int argc, char *argv[])
{
  sim::Common c;
  sim::Simulation s(c);

  set_parameters(s.common);
  sim::process_command_line(c, argc, argv);
  s.common.set_defaults_not_yet_set();
  s.common.convert_probabilities_to_time_period("RISK_DEATH_MALE", sim::YEAR);
  s.common.convert_probabilities_to_time_period("RISK_DEATH_FEMALE", sim::YEAR);
  s.common.convert_probabilities_to_time_period("HIV_RISK_DEATH",
						sim::YEAR, 1, 2);
  for (unsigned i = 0; i < c("NUM_SIMULATIONS"); ++i) {
    sim::Common common = s.common;
    sim::Simulation simulation(common);
    simulation.init({sim::advanceTimeEvent, cd4Event, deathEvent});
    simulation.simulate();
    report(simulation, i);
  }

  return 0;
}
