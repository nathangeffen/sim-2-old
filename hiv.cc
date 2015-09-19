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

 Deterministic population growth E1
 Deterministic population deaths by HIV and ARV status E2
     (5-stage Granich model)
 Deterministic disease progression E3
 Stochastic HIV infection (including by risk strata) 5 E6
 Deterministic ARVs E5

 ### Model C

 Deterministic population growth E1
 Stochastic population deaths by HIV and ARV status (infection time) E7
 Deterministic disease progression E3
 Stochastic HIV infection (including by risk strata) 5 E6
 Deterministic ARVs E5

 ### Model D

 Deterministic population growth E1
 Stochastic population deaths by HIV and ARV status (infection time) E7
 Deterministic disease progression E3
 Stochastic HIV infection (including by risk strata) 5 E6
 Stochastic ARVs E8

 ### Model E

 Stochastic population growth E9
 Stochastic population deaths by HIV and ARV status (infection time) E7
 Stochastic disease progression E10
 Stochastic HIV infection (including by risk strata) 5 E6
 Stochastic ARVs E8

 ### Model F

 Stochastic population immigration
 Stochastic population emigration
 Stochastic population births
 Stochastic population deaths by HIV and ARV status (infection time)
 Stochastic HIV infection using Leigh's method E11

 ### Model G

 Stochastic population immigration
 Stochastic population emigration
 Stochastic population births
 Stochastic population deaths by HIV and ARV status (infection time)
 Stochastic random matching for HIV infection E12


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

class HIVAgent : public Agent {
public:
  HIVAgent(Context &c) : Agent(c) {};
  int hiv;
  int risk;
  double hiv_infection_date;
  double on_arvs_date = 0.0;
  double cd4;
  bool circumcised;
  double orientation;
  double riskiness;
  double expected_encounters;
  unsigned max_iteration_encounters;
  unsigned iteration_encounters;
  unsigned total_encounters;
  double cluster_val;
  std::vector<HIVAgent *> partners;
  inline double cluster()
  {
    return riskiness;
  }
  friend std::ostream& operator<<(std::ostream& os, const HIVAgent& agent)
  {
    os << " ID: " << agent.id << " Birth: " << agent.dob
       << " Death: " << agent.dod
       << " HIV: " << agent.hiv << " HIV infection date: "
       << agent.hiv_infection_date << " CD4: " << agent.cd4;
    return os;
  }
};

HIVAgent *create_hiv_agent(Context &c)
{
  std::uniform_real_distribution<double> uni;
  beta_distribution<> beta(c("BETA_DIST_ALPHA"), c("BETA_DIST_BETA"));
  HIVAgent *a = new HIVAgent(c);
  a->sex = uni(rng) < c("PROB_MALE") ?
		      MALE : FEMALE;
  { // dob
    std::uniform_real_distribution<double>
      uni_age(c("EARLIEST_BIRTH_DATE"),
	      c("LATEST_BIRTH_DATE"));
    a->dob = uni_age(rng);
  }

  a->cd4 = 1000;
  a->hiv = 0;
  if (uni(rng) < c("HIV_PREVALENCE")) {
    a->hiv = 1;
    std::uniform_real_distribution<double>
      uni_age(c("EARLIEST_BIRTH_DATE"),
	      c("LATEST_BIRTH_DATE"));
  }
  a->riskiness = beta(rng);
  a->risk = uni(rng) < c("HIGH_RISK_PROPORTION") ? 1 : 0;
  a->expected_encounters = a->risk == 0 ?
    c("NUM_PARTNERSHIPS_LOW") : c("NUM_PARTNERSHIPS_HIGH");
  a->total_encounters = 0;
  a->orientation = 1.0;


  if (uni(rng) < c("PROB_CIRCUMCISED"))
    a->circumcised = true;

  return a;
}

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

void setDefaultParameters(sim::Simulation &s)
{
  s.context.set("NUM_AGENTS", {1000.0});
  s.context.set("INITIAL_AGE", {15.0});
  s.context.set("BETA_DIST_ALPHA", { 2.0 });
  s.context.set("BETA_DIST_BETA", { 2.0 });
}

/* E1 */

void increasePopulationDeterministic(sim::Simulation &s)
{
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

/*
 * E6 - Infection simple stochastic
 *
 */


void infectionRiskSimpleStochastic(sim::Simulation &s)
{
  size_t total_infections = 0;
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
  // Low risk
  double num_partnerships_low = s.context("NUM_PARTNERSHIPS_LOW");
  double risk_not_infected_low = pow(risk_not_infected_individual_act,
				     num_partnerships_low);
  double risk_infected_low = 1.0 - risk_not_infected_low;
  risk_infected_low = sim::time_correct_prob(risk_infected_low,
					     time_period,
					     s.time_step);


  for (auto a = s.agents.begin(); a != s.agents.end(); ++a) {
    HIVAgent *agent = (HIVAgent *) *a;
    if (agent->risk == 1 && agent->hiv == 0 && is_event(risk_infected_high)) {
      agent->hiv = 1;
      agent->hiv_infection_date = s.current_date;
      ++total_infections;
    } else if (agent->risk == 0 && agent->hiv == 0 &&
	       is_event(risk_infected_low)) {
      agent->hiv = 1;
      agent->hiv_infection_date = s.current_date;
      ++total_infections;
    }
  }

  s.context.tracking.new_infections.push_back(total_infections);
}


/*
 * E7
 *
 * Stage 0: HIV-negative
 * Stage 1 to 4: WHO stages
 * Stage 5: ARVs
 */

void decreasePopulationStochastic(sim::Simulation &s)
{
  double decline_rates[6];
  for (int stage = 0; stage < 6; ++stage)
    decline_rates[stage] = s.context("DECLINE_RATE", stage);

  for (auto & a : s.agents) {
    HIVAgent *agent = (HIVAgent *) a;
    if (is_event(decline_rates[agent->hiv]))
      agent->die(s, "DETERMINISTIC");
  }

  s.remove_dead_agents();
}

/*
 * E8
 *
 * People in each HIV stage move to ARVs with differing probabilities.  People
 * on ARVs can also stop treatment or have viral rebound and move back to stage
 * 3.
 *
 * ARV_RATE 0, 1, 2 and 3 are the rates of moving onto ARVs for WHO stages 1, 2,
 * 3 and 4 respectively. ARV_RATE 5 is the rate for moving off ARVs.
 */

void arvsStochastic(sim::Simulation &s)
{
  double arv_rates[6];
  for (int stage = 1; stage < 6; ++stage)
    arv_rates[stage] = s.context("ARV_RATE", stage - 1);

  for (auto & a : s.agents) {
    HIVAgent *agent = (HIVAgent *) a;
    if (agent->hiv > 0 && is_event(arv_rates[agent->hiv])) {
      if (agent->hiv < 5)
	agent->hiv = 5;
      else
	agent->hiv = 3;
    }
  }
}

/* E9 */

void increasePopulationStochastic(sim::Simulation &s)
{
  double growth_rate = s.context("GROWTH");
  size_t num_agents = s.context("_AGENTS");
  double initial_age = s.context("INITIAL_AGE");

  for (size_t i = 0; i < num_agents; ++i) {
    if (is_event(growth_rate)) {
      HIVAgent *a = create_hiv_agent(s.context);
      a->dob = s.current_date - initial_age;
      a->hiv = 0;
      s.agents.push_back(a);
    }
  }
}


/*
 * E10
 *
 * People move from HIV stage 1 through to HIV infection WHO stage 4,
 *
 */

void diseaseProgressionStochastic(sim::Simulation &s)
{
  double change_rates[5];
  for (int stage = 1; stage < 4; ++stage)
    change_rates[stage] = s.context("DISEASE_CHANGE_RATE", stage - 1);

  for (auto & a : s.agents) {
    HIVAgent *agent = (HIVAgent *) a;
    if (agent->hiv > 0 && agent->hiv < 4 && is_event(change_rates[agent->hiv]))
      ++agent->hiv;
  }
}

/* Support functions for partner matching algorithms. */

unsigned expected_encounters(const Simulation &s, const HIVAgent *a)
{
  return a->expected_encounters;
}

unsigned poisson_encounters(const Simulation &s, const HIVAgent *a)
{
  std::poisson_distribution<unsigned> dist(a->expected_encounters);
  return dist(sim::rng);
}

HIVAgent* use_existing_partners(HIVAgent *agent)
{
  HIVAgent *partner;
  do {
    std::uniform_int_distribution<unsigned>
      dist(0, agent->partners.size() - 1);
    unsigned index = dist(sim::rng);
    partner = (HIVAgent *) agent->partners[index];
    if (partner->alive == false) {
      partner = NULL;
      agent->partners[index] = agent->partners.back();
      agent->partners.pop_back();
    }
  } while(partner == NULL && agent->partners.size() > 0);
  return partner;
}


HIVAgent *random_partner(const Simulation &s,
			 HIVAgent *a,
			 std::vector<Agent *>::const_iterator from,
			 std::vector<Agent *>::const_iterator to)
{
  std::uniform_int_distribution<unsigned> dist(0, to - from - 1);
  return (HIVAgent *) from[dist(sim::rng)];
}

HIVAgent *random_sticky_partner(const Simulation &s,
				HIVAgent *a,
				std::vector<Agent *>::const_iterator from,
				std::vector<Agent *>::const_iterator to)
{
  HIVAgent *partner = NULL;
  bool change_partners = is_event(a->riskiness);
  // Try to match. "attempts" ensures no endless loop (unlikely)
  unsigned attempts = 0;
  while (partner == NULL && attempts < 5) { // 5 is arbitrary choice
    if (a->partners.size() == 0 || change_partners) {
      // Find a new partner
      std::uniform_int_distribution<unsigned> dist(0, to - from - 1);
      auto p = from[dist(sim::rng)];
      auto index = std::find(a->partners.begin(), a->partners.end(), p);
      if (index == a->partners.end()) { // Found a new partner
	partner = (HIVAgent *) p;
	a->partners.push_back(partner);
	partner->partners.push_back(a);
      }
    } else {
      partner = use_existing_partners(a);
    }
    ++attempts;
  }
  return partner;
}

double distance_simple(const Simulation &s,
		       const HIVAgent *a,
		       const HIVAgent *b)
{
  double distance = 0.0;
  if (a->sex == b->sex)
    distance += 50000.0;
  distance += fabs(a->riskiness - b->riskiness);
  return distance;
}

class SampleKPartners {
public:
  SampleKPartners(unsigned k,
		  std::function<double(const Simulation &,
				       const HIVAgent *,
				       const HIVAgent *)> distance_func)
    : k_(k), distance_func_(distance_func) {}
  HIVAgent * operator()(const Simulation &s,
			HIVAgent *a,
			std::vector<Agent *>::const_iterator from,
			std::vector<Agent *>::const_iterator to)
  {
    HIVAgent *partner = NULL;
    bool change_partners = is_event(a->riskiness);
    // Try to match. *attempts* ensures no endless loop (unlikely)
    unsigned attempts = 0;
    while (partner == NULL && attempts < 5) { // 5 is arbitrary choice
      if (a->partners.size() == 0 || change_partners) {
	// Find a new partner
	HIVAgent *best_partner = NULL;
	double best_distance = std::numeric_limits<double>::max();
	unsigned i = 0;
	auto p = to - 1;
	for (; p >=from && i < k_; --p, ++i) {
	  double d = distance_func_(s, a, (HIVAgent *) *p);
	  if (d < best_distance) {
	    auto index = std::find(a->partners.begin(), a->partners.end(), *p);
	    if (index == a->partners.end()) {
	      best_distance = d;
	      best_partner = (HIVAgent *) *p;
	    }
	  }
	}
	partner = best_partner;
      } else {
	partner = use_existing_partners(a);
      }
      ++attempts;
    }
    return partner;
  }
private:
  unsigned k_;
  std::function<double(const Simulation &,
		       const HIVAgent *, const HIVAgent *)> distance_func_;
};



void weighted_shuffle(sim::Simulation & s)
{
  std::uniform_real_distribution<double> uni;
  for (auto & a : s.agents) {
    HIVAgent *agent = (HIVAgent *) a;
    agent->cluster_val = agent->cluster() * uni(rng);
  }
  std::sort(s.agents.begin(), s.agents.end(),
	    [](const Agent *a, const Agent *b) {
	      HIVAgent *a1 = (HIVAgent *) a;
	      HIVAgent *b1 = (HIVAgent *) b;
	      return a1->cluster_val > b1->cluster_val;
	    });
}

void cluster_shuffle(sim::Simulation &s)
{
  size_t clusters = s.context("CLUSTER_SIZE");
  size_t cluster_size = s.agents.size() / clusters;
  std::sort(s.agents.begin(), s.agents.end(),
	    [](const Agent *a, const Agent *b) {
	      HIVAgent *a1 = (HIVAgent *) a;
	      HIVAgent *b1 = (HIVAgent *) b;
	      return a1->cluster() > b1->cluster();
	    });
  for (size_t c = 0; c < clusters; ++c) {
    auto first = c * cluster_size;
    auto last = first + cluster_size;
    shuffle(s.agents.begin() + first, s.agents.begin() + last, sim::rng);
  }
}


/*
 * E12
 *
 * Matches individuals randomly for sex.
 *
 * HOW THIS WORKS
 *
 * Shuffle agents for stochastic effects
 * Determine how many times each individual will have sex
 * Iterate backwards through agents
 *   For each agent, repeat for number of times s/he has sex
 *      Determine if partner change and if so, randomly find partner
 *       (being careful not to select one that's already selected)
 *      If no partner change, randomly select existing partner
 *       (being careful not to mate with a dead one)
 *
 * Many shortcomings, but particularly that although
 * partners increase, unless they die, they never decrease.
 */


class InfectionRiskRandomMatching
{
public:
  InfectionRiskRandomMatching(std::function<unsigned(const Simulation &,
						     const HIVAgent *)>
			      get_num_partners_func = expected_encounters,
			      std::function<HIVAgent *
			      (const Simulation &s,
			       HIVAgent *a,
			       std::vector<Agent *>::const_iterator from,
			       std::vector<Agent *>::const_iterator to)>
			      get_partner_func = random_partner,
			      std::function<void(sim::Simulation &)>
			      shuffle_func = [](sim::Simulation &s) {
				std::shuffle(s.agents.begin(), s.agents.end(),
					     sim::rng);
			      },
			      bool reverse_agents = false)
  {
    get_num_partners_func_ = get_num_partners_func;
    get_partner_func_ = get_partner_func;
    shuffle_func_ = shuffle_func;
    reverse_agents_ = reverse_agents;
  }

  void operator()(sim::Simulation &s)
  {
    shuffle_func_(s);
    if (reverse_agents_)
      reverse(s.agents.begin(), s.agents.end());

    // Determine how many times each individual will have sex
    unsigned total = 0, encounters = 0, infections = 0;
    for (size_t i = 0; i < s.agents.size(); ++i) {
      HIVAgent *a = (HIVAgent *) s.agents[i];
      a->max_iteration_encounters = get_num_partners_func_(s, a);
      a->iteration_encounters = 0;
      total += a->max_iteration_encounters;
    }

    for (size_t i = s.agents.size() - 1; i > 0 && encounters < total / 2; --i) {
      HIVAgent *a = (HIVAgent *) s.agents[i];
      for (;a->iteration_encounters < a->max_iteration_encounters;
	   ++a->iteration_encounters) {
	HIVAgent *partner = get_partner_func_(s, a,
					      s.agents.begin(),
					      s.agents.begin() + i);
	if (partner) { // unlikely that this is false
	  ++a->total_encounters;
	  ++partner->total_encounters;
	  ++partner->iteration_encounters;
	  ++encounters;
	  // Check for discordance
	  // a infected, partner not
	  if (a->hiv > 0 && a->hiv < 5 && partner->hiv == 0) {
	    // infection event
	    if (is_event(s.context("RISK_TRANSMISSION_ACT"))) {
	      partner->hiv = 1;
	      partner->hiv_infection_date = s.current_date;
	      ++infections;
	    }
	    // partner infected, a not
	  } else if (partner->hiv > 0 && partner->hiv < 5 && a->hiv == 0) {
	    // infection event
	    if (is_event(s.context("RISK_TRANSMISSION_ACT"))) {
	      a->hiv = 1;
	      a->hiv_infection_date = s.current_date;
	      ++infections;
	    }
	  }
	} else {
	  std::cout << "WARNING: Unlikely non-match error occurred with agent "
		    << a->id << std::endl;
	}
      }
    }
    s.context.tracking.new_infections.push_back(infections);
  }
private:
  std::function<unsigned(const Simulation &,
			 const HIVAgent *)> get_num_partners_func_;
  std::function<HIVAgent * (const Simulation &,
			    HIVAgent *,
			    std::vector<Agent *>::const_iterator,
			    std::vector<Agent *>::const_iterator)>
  get_partner_func_;
  std::function<void(sim::Simulation &)> shuffle_func_;
  bool reverse_agents_;
};



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
    HIVAgent *agent = new HIVAgent(s.context);
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
		  .beforeAllSimulations(setDefaultParameters)
		  .afterEachSimulation([&t](sim::Simulation &s) {
		      TESTEQ(t, s.agents.size(), 1220, "Population growth");
		    })
		  .agentCreate(create_hiv_agent)
		  .timeAdjust("GROWTH", 1.0, 0, 1, sim::PROBABILITY)
		  .parameter("HIGH_RISK_PROPORTION", {0.5 } )
		  .parameter("TIME_STEP", {1.0} )
		  .parameter("GROWTH", {0.01} )
		  .parameter("RISK_SUSCEPTIBLE", {0.0} )
		  .parameter("RISK_SUSCEPTIBLE_STDEV", {0.0} )
		  .parameter("RISK_UNINFECTIOUS", {0.0} )
		  .parameter("RISK_UNINFECTIOUS_STDEV", {0.0} )
		  .parameter("PARTNER_FORMATION_RATE", {0.0} )
		  .parameter("PARTNER_FORMATION_RATE_STDEV", {0.0} )
		  .parameter("NUM_PARTNERSHIPS_LOW", {0} )
		  .parameter("NUM_PARTNERSHIPS_HIGH", {6} )
		  .parameter("HIGH_RISK_PROPORTION", {0.5 } )
		  .parameter("PREVALENCE", {0.0} )).simulate();

  std::cout << "Testing E1 time step half year" << std::endl;
  sim::Simulation(sim::Options()
		  .events({sim::advanceTimeEvent
			, calcVariablesEvent
			, increasePopulationDeterministic})
		  .beforeAllSimulations(setDefaultParameters)
		  .afterEachSimulation([&t](sim::Simulation &s) {
		      TESTEQ(t, s.agents.size(), 1221, "Population growth");
		    })
		  .agentCreate(create_hiv_agent)
		  .timeAdjust("GROWTH", 1.0, 0, 1, sim::PROBABILITY)
		  .parameter("HIGH_RISK_PROPORTION", {0.5 } )
		  .parameter("TIME_STEP", {0.5} )
		  .parameter("GROWTH", {0.01} )
		  .parameter("RISK_SUSCEPTIBLE", {0.0} )
		  .parameter("RISK_SUSCEPTIBLE_STDEV", {0.0} )
		  .parameter("RISK_UNINFECTIOUS", {0.0} )
		  .parameter("RISK_UNINFECTIOUS_STDEV", {0.0} )
		  .parameter("PARTNER_FORMATION_RATE", {0.0} )
		  .parameter("PARTNER_FORMATION_RATE_STDEV", {0.0} )
		  .parameter("NUM_PARTNERSHIPS_LOW", {0} )
		  .parameter("NUM_PARTNERSHIPS_HIGH", {6} )
		  .parameter("HIGH_RISK_PROPORTION", {0.5 } )
		  .parameter("PREVALENCE", {0.0} )).simulate();

  std::cout << "Testing E1 time step 1 week" << std::endl;
  sim::Simulation(sim::Options()
		  .events({sim::advanceTimeEvent
			, calcVariablesEvent
			, increasePopulationDeterministic})
		  .beforeAllSimulations(setDefaultParameters)
		  .afterEachSimulation([&t](sim::Simulation &s) {
		      TESTEQ(t, s.agents.size(), 1222, "Population growth");
		    })
		  .agentCreate(create_hiv_agent)
		  .timeAdjust("GROWTH", 1.0, 0, 1, sim::PROBABILITY)
		  .parameter("NUM_AGENTS", {1000.0} )
		  .parameter("HIGH_RISK_PROPORTION", {0.5 } )
		  .parameter("TIME_STEP", {sim:WEEK} )
		  .parameter("GROWTH", {0.01} )
		  .parameter("RISK_SUSCEPTIBLE", {0.0} )
		  .parameter("RISK_SUSCEPTIBLE_STDEV", {0.0} )
		  .parameter("RISK_UNINFECTIOUS", {0.0} )
		  .parameter("RISK_UNINFECTIOUS_STDEV", {0.0} )
		  .parameter("PARTNER_FORMATION_RATE", {0.0} )
		  .parameter("PARTNER_FORMATION_RATE_STDEV", {0.0} )
		  .parameter("NUM_PARTNERSHIPS_LOW", {0} )
		  .parameter("NUM_PARTNERSHIPS_HIGH", {6} )
		  .parameter("HIGH_RISK_PROPORTION", {0.5 } )
		  .parameter("BETA_DIST_ALPHA", { 2.0 })
		  .parameter("BETA_DIST_BETA", { 2.0 })
		  .parameter("PREVALENCE", {0.0} )).simulate();


  // Test E2 time step 1
  std::cout << "Testing E2 time step 1" << std::endl;
  sim::Simulation(sim::Options()
		  .events({sim::advanceTimeEvent
			, calcVariablesEvent
			, decreasePopulationDeterministic})
		  .beforeAllSimulations(setDefaultParameters)
		  .afterEachSimulation([&t](sim::Simulation &s) {
		      TESTEQ(t, s.agents.size(), 817, "Population decline");
		    })
		  .agentCreate(create_hiv_agent)
		  .timeAdjust("DECLINE_RATE", 1.0, 0, 1, sim::PROBABILITY)
		  .parameter("HIGH_RISK_PROPORTION", {0.5 } )
		  .parameter("DECLINE_RATE",
			     {0.01, 0.02, 0.03, 0.1, 0.5, 0.015} )
		  .parameter("TIME_STEP", {1.0} )
		  .parameter("HIV_PREVALENCE", {0.0, 0.0, 0.0, 0.0} )
		  .parameter("NUM_PARTNERSHIPS_LOW", {0} )
		  .parameter("NUM_PARTNERSHIPS_HIGH", {6} )
		  .parameter("HIGH_RISK_PROPORTION", {0.5 } )).simulate();

  // Test E2 time step half-year
  std::cout << "Testing E2 time step half year" << std::endl;
  sim::Simulation(sim::Options()
		  .events({sim::advanceTimeEvent
			, calcVariablesEvent
			, decreasePopulationDeterministic})
		  .beforeAllSimulations(setDefaultParameters)
		  .afterEachSimulation([&t](sim::Simulation &s) {
		      TESTEQ(t, s.agents.size(), 817, "Population decline");
		    })
		  .agentCreate(create_hiv_agent)
		  .timeAdjust("DECLINE_RATE", 1.0, 0, 1, sim::PROBABILITY)
		  .parameter("HIGH_RISK_PROPORTION", {0.5 } )
		  .parameter("DECLINE_RATE",
			     {0.01, 0.02, 0.03, 0.1, 0.5, 0.015} )
		  .parameter("TIME_STEP", { 0.5 } )
		  .parameter("HIV_PREVALENCE", {0.0, 0.0, 0.0, 0.0} )
		  .parameter("NUM_PARTNERSHIPS_LOW", {0} )
		  .parameter("NUM_PARTNERSHIPS_HIGH", {6} )
		  .parameter("HIGH_RISK_PROPORTION", {0.5 } )).simulate();

  // Test E2 time step week
  std::cout << "Testing E2 time step week" << std::endl;
  sim::Simulation(sim::Options()
		  .events({sim::advanceTimeEvent
			, calcVariablesEvent
			, decreasePopulationDeterministic})
		  .beforeAllSimulations(setDefaultParameters)
		  .afterEachSimulation([&t](sim::Simulation &s) {
		      TESTEQ(t, s.agents.size(), 817, "Population decline");
		    })
		  .agentCreate(create_hiv_agent)
		  .timeAdjust("DECLINE_RATE", 1.0, 0, 1, sim::PROBABILITY)
		  .parameter("HIGH_RISK_PROPORTION", {0.5 } )
		  .parameter("DECLINE_RATE",
			     {0.01, 0.02, 0.03, 0.1, 0.5, 0.015} )
		  .parameter("TIME_STEP", { sim::WEEK } )
		  .parameter("HIV_PREVALENCE", {0.0, 0.0, 0.0, 0.0} )
		  .parameter("NUM_PARTNERSHIPS_LOW", {0} )
		  .parameter("NUM_PARTNERSHIPS_HIGH", {6} )
		  .parameter("HIGH_RISK_PROPORTION", {0.5 } )).simulate();


  // Test E3
  std::cout << "Testing E3 time step year" << std::endl;
  sim::Simulation(sim::Options()
		  .events({sim::advanceTimeEvent
			, calcVariablesEvent
			, diseaseProgressionDeterministic})
		  .beforeAllSimulations(setDefaultParameters)
		  .afterEachSimulation([&t](sim::Simulation &s) {
		      unsigned stage_1 =
			count_if(s.agents.begin(), s.agents.end(),
				 [](const sim::Agent *a) {
				   const HIVAgent *b = (HIVAgent *) a;
				   return b->hiv == 1;
				 });
		      TESTEQ(t, stage_1, 121, "Infection change");
		    })
		  .agentCreate(create_hiv_agent)
		  .timeAdjust("DISEASE_CHANGE_RATE", 1.0, 0, 1,
			      sim::PROBABILITY)
		  .parameter("HIGH_RISK_PROPORTION", {0.5 } )
		  .parameter("DISEASE_CHANGE_RATE", {0.1, 0.2, 0.3} )
		  .parameter("TIME_STEP", {1.0 } )
		  .parameter("HIV_PREVALENCE", {1.0, 0.0, 0.0, 0.0} )
		  .parameter("NUM_PARTNERSHIPS_LOW", {0} )
		  .parameter("NUM_PARTNERSHIPS_HIGH", {6} )
		  .parameter("HIGH_RISK_PROPORTION", {0.5 } )).simulate();

  std::cout << "Testing E3 time step half-year" << std::endl;
  sim::Simulation(sim::Options()
		  .events({sim::advanceTimeEvent
			, calcVariablesEvent
			, diseaseProgressionDeterministic})
		  .beforeAllSimulations(setDefaultParameters)
		  .afterEachSimulation([&t](sim::Simulation &s) {
		      unsigned stage_1 =
			count_if(s.agents.begin(), s.agents.end(),
				 [](const sim::Agent *a) {
				   const HIVAgent *b = (HIVAgent *) a;
				   return b->hiv == 1;
				 });
		      TESTEQ(t, stage_1, 121, "Infection change");
		    })
		  .agentCreate(create_hiv_agent)
		  .timeAdjust("DISEASE_CHANGE_RATE", 1.0, 0, 1,
			      sim::PROBABILITY)
		  .parameter("HIGH_RISK_PROPORTION", {0.5 } )
		  .parameter("DISEASE_CHANGE_RATE", {0.1, 0.2, 0.3} )
		  .parameter("TIME_STEP", {0.5 } )
		  .parameter("HIV_PREVALENCE", {1.0, 0.0, 0.0, 0.0} )
		  .parameter("NUM_PARTNERSHIPS_LOW", {0} )
		  .parameter("NUM_PARTNERSHIPS_HIGH", {6} )
		  .parameter("HIGH_RISK_PROPORTION", {0.5 } )).simulate();

  std::cout << "Testing E3 time step week." << std::endl;
  sim::Simulation(sim::Options()
		  .events({sim::advanceTimeEvent
			, calcVariablesEvent
			, diseaseProgressionDeterministic})
		  .beforeAllSimulations(setDefaultParameters)
		  .afterEachSimulation([&t](sim::Simulation &s) {
		      unsigned stage_1 =
			count_if(s.agents.begin(), s.agents.end(),
				 [](const sim::Agent *a) {
				   const HIVAgent *b = (HIVAgent *) a;
				   return b->hiv == 1;
				 });
		      TESTEQ(t, stage_1, 121, "Infection change");
		    })
		  .agentCreate(create_hiv_agent)
		  .timeAdjust("DISEASE_CHANGE_RATE", 1.0, 0, 1,
			      sim::PROBABILITY)
		  .parameter("HIGH_RISK_PROPORTION", {0.5 } )
		  .parameter("DISEASE_CHANGE_RATE", {0.1, 0.2, 0.3} )
		  .parameter("TIME_STEP", {sim::WEEK} )
		  .parameter("HIV_PREVALENCE", {1.0, 0.0, 0.0, 0.0} )
		  .parameter("NUM_PARTNERSHIPS_LOW", {0} )
		  .parameter("NUM_PARTNERSHIPS_HIGH", {6} )
		  .parameter("HIGH_RISK_PROPORTION", {0.5 } )).simulate();



  // Test E4 - I can't figure out how to time correct these parameters. The
  // equations are too hard. See calculations.ods.

  std::cout << "Testing E4 one year with none on ARVs" << std::endl;
  sim::Simulation(sim::Options()
		  .events({sim::advanceTimeEvent
			, calcVariablesEvent
			, infectionRiskDeterministic})
		  .beforeAllSimulations(setDefaultParameters)
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
		  .parameter("HIGH_RISK_PROPORTION", {0.5 } )
		  .parameter("TIME_STEP", {1.0} )
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
		  .beforeAllSimulations(setDefaultParameters)
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
		  .parameter("HIGH_RISK_PROPORTION", {0.5 } )
		  .parameter("TIME_STEP", {1.0} )
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
		  .beforeAllSimulations(setDefaultParameters)
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
		  .parameter("HIGH_RISK_PROPORTION", {0.5 } )
		  .parameter("TIME_STEP", {1.0} )
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
		  .beforeAllSimulations(setDefaultParameters)
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
  		  .parameter("TIME_STEP", {1.0 } )
  		  .parameter("HIV_PREVALENCE", {0.0, 0.0, 0.0, 0.0} )
  		  .parameter("HIV_PREVALENCE_ARVS", {0.5} )
  		  .parameter("HIGH_RISK_PROPORTION", {0.5 } )
		  .parameter("NUM_PARTNERSHIPS_LOW", {0} )
		  .parameter("NUM_PARTNERSHIPS_HIGH", {6} )).simulate();


  std::cout << "Testing E5 week - ARVS to WHO stage only" << std::endl;
  sim::Simulation(sim::Options()
  		  .events({sim::advanceTimeEvent
  			, calcVariablesEvent
  			, arvsDeterministic})
		  .beforeAllSimulations(setDefaultParameters)
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
  		  .parameter("TIME_STEP", { sim::WEEK } )
  		  .parameter("HIV_PREVALENCE", {0.0, 0.0, 0.0, 0.0} )
  		  .parameter("HIV_PREVALENCE_ARVS", {0.5} )
  		  .parameter("HIGH_RISK_PROPORTION", {0.5 } )
		  .parameter("NUM_PARTNERSHIPS_LOW", {0} )
		  .parameter("NUM_PARTNERSHIPS_HIGH", {6} )).simulate();

  std::cout << "Testing E5 year - WHO stage only to ARVs" << std::endl;
  sim::Simulation(sim::Options()
  		  .events({sim::advanceTimeEvent
  			, calcVariablesEvent
  			, arvsDeterministic})
		  .beforeAllSimulations(setDefaultParameters)
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
  		  .parameter("TIME_STEP", { sim::YEAR } )
  		  .parameter("HIV_PREVALENCE", {0.1, 0.2, 0.3, 0.1} )
  		  .parameter("HIV_PREVALENCE_ARVS", {0.1} )
  		  .parameter("HIGH_RISK_PROPORTION", {0.5 } )
		  .parameter("NUM_PARTNERSHIPS_LOW", {0} )
		  .parameter("NUM_PARTNERSHIPS_HIGH", {6} )).simulate();

  std::cout << "Testing E6 one year with none on ARVs" << std::endl;
  const size_t simulations = 100;
  std::vector<double> measures(simulations, 0.0);
  sim::Simulation(sim::Options()
		  .events({sim::advanceTimeEvent
			, calcVariablesEvent
			, infectionRiskSimpleStochastic})
		  .beforeAllSimulations(setDefaultParameters)
		  .afterEachSimulation([&measures](sim::Simulation &s) {
		      unsigned susceptible =
			count_if(s.agents.begin(), s.agents.end(),
				 [](const sim::Agent *a) {
				   const HIVAgent *b = (HIVAgent *) a;
				   return b->hiv == 0;
				 });
		      measures[s.simulation_num - 1] = (double) susceptible;
		    })
		  .allAgentsCreate(createDeterministicAgents)
		  .parameter("HIGH_RISK_PROPORTION", {0.5 } )
		  .parameter("TIME_STEP", {1.0} )
		  .parameter("HIGH_RISK_PROPORTION", {0.5 } )
		  .parameter("INFECTION_PARAMETERS_TIME", {1.0})
		  .parameter("FORCE_INFECTION", {0.05} )
		  .parameter("NUM_PARTNERSHIPS_LOW", {0} )
		  .parameter("NUM_PARTNERSHIPS_HIGH", {6} )
		  .parameter("HIV_PREVALENCE", {0.1, 0.00, 0.00, 0.00} )
		  .parameter("NUM_SIMULATIONS", {simulations})
		  .parameter("HIV_PREVALENCE_ARVS", {0.00} ))
		  .simulate();
  TESTRANGE(t, sim::mean(measures), 583.0, 595.0, "Mean of stochastic run.");

  measures.clear();
  measures.resize(simulations);
  std::cout << "Testing E6 one year with some on ARVs" << std::endl;
  sim::Simulation(sim::Options()
		  .events({sim::advanceTimeEvent
			, calcVariablesEvent
			, infectionRiskSimpleStochastic})
		  .beforeAllSimulations(setDefaultParameters)
		  .afterEachSimulation([&measures](sim::Simulation &s) {
		      unsigned susceptible =
			count_if(s.agents.begin(), s.agents.end(),
				 [](const sim::Agent *a) {
				   const HIVAgent *b = (HIVAgent *) a;
				   return b->hiv == 0;
				 });
		      measures[s.simulation_num - 1] = (double) susceptible;
		    })
		  .allAgentsCreate(createDeterministicAgents)
		  .parameter("HIGH_RISK_PROPORTION", {0.5 } )
		  .parameter("TIME_STEP", {1.0} )
		  .parameter("HIGH_RISK_PROPORTION", {0.5 } )
		  .parameter("INFECTION_PARAMETERS_TIME", {1.0})
		  .parameter("FORCE_INFECTION", {0.05} )
		  .parameter("NUM_PARTNERSHIPS_LOW", {0} )
		  .parameter("NUM_PARTNERSHIPS_HIGH", {6} )
		  .parameter("HIV_PREVALENCE", {0.05, 0.00, 0.00, 0.00} )
		  .parameter("NUM_SIMULATIONS", {simulations})
		  .parameter("HIV_PREVALENCE_ARVS", {0.1} ))
		  .simulate();
  TESTRANGE(t, sim::mean(measures), 661.0, 669.0, "Mean of stochastic run.");

  measures.clear();
  measures.resize(simulations);
  std::cout << "Testing E7 one year" << std::endl;
  sim::Simulation(sim::Options()
  		  .events({sim::advanceTimeEvent
  			, calcVariablesEvent
  			, decreasePopulationStochastic})
		  .beforeAllSimulations(setDefaultParameters)
  		  .afterEachSimulation([&measures](sim::Simulation &s) {
		      measures[s.simulation_num - 1] = (double) s.agents.size();
  		    })
  		  .agentCreate(create_hiv_agent)
  		  .timeAdjust("DECLINE_RATE", 1.0, 0, 1, sim::PROBABILITY)
  		  .parameter("HIGH_RISK_PROPORTION", {0.5 } )
  		  .parameter("DECLINE_RATE",
  			     {0.01, 0.02, 0.03, 0.1, 0.5, 0.015} )
  		  .parameter("TIME_STEP", {1.0} )
  		  .parameter("HIV_PREVALENCE", {0.0, 0.0, 0.0, 0.0} )
  		  .parameter("NUM_SIMULATIONS", {simulations})
		  .parameter("NUM_PARTNERSHIPS_LOW", {0} )
		  .parameter("NUM_PARTNERSHIPS_HIGH", {6} )
  		  .parameter("HIGH_RISK_PROPORTION", {0.5 } )).simulate();
  TESTRANGE(t, sim::mean(measures), 810.0, 822.0, "Mean of stochastic run.");

  std::cout << "Testing E8 one year - ARVS to WHO stage only" << std::endl;
  measures.clear();
  measures.resize(simulations);
  sim::Simulation(sim::Options()
  		  .events({sim::advanceTimeEvent
  			, calcVariablesEvent
  			, arvsStochastic})
		  .beforeAllSimulations(setDefaultParameters)
  		  .afterEachSimulation([&measures](sim::Simulation &s) {
  		      unsigned on_arvs =
  			count_if(s.agents.begin(), s.agents.end(),
  				 [](const sim::Agent *a) {
  				   const HIVAgent *b = (HIVAgent *) a;
  				   return b->hiv == 5;
  				 });
		      measures[s.simulation_num - 1] = on_arvs;
  		    })
  		  .allAgentsCreate(createDeterministicAgents)
  		  .timeAdjust("ARV_RATE", 1.0, 0, 1, sim::PROBABILITY)
  		  .parameter("ARV_RATE", {0.0, 0.0, 0.0, 0.0, 0.1} )
  		  .parameter("TIME_STEP", {1.0 } )
  		  .parameter("HIV_PREVALENCE", {0.0, 0.0, 0.0, 0.0} )
  		  .parameter("HIV_PREVALENCE_ARVS", {0.5} )
  		  .parameter("HIGH_RISK_PROPORTION", {0.5 } )
  		  .parameter("NUM_SIMULATIONS", {simulations})
		  .parameter("NUM_PARTNERSHIPS_LOW", {0} )
		  .parameter("NUM_PARTNERSHIPS_HIGH", {6} )).simulate();
  TESTRANGE(t, sim::mean(measures), 55, 65, "Mean of stochastic run.");

  std::cout << "Testing E8 half-year - ARVs to WHO stage only" << std::endl;
  measures.clear();
  measures.resize(simulations);
  sim::Simulation(sim::Options()
  		  .events({sim::advanceTimeEvent
  			, calcVariablesEvent
  			, arvsStochastic})
		  .beforeAllSimulations(setDefaultParameters)
  		  .afterEachSimulation([&measures](sim::Simulation &s) {
  		      unsigned on_arvs =
  			count_if(s.agents.begin(), s.agents.end(),
  				 [](const sim::Agent *a) {
  				   const HIVAgent *b = (HIVAgent *) a;
  				   return b->hiv == 5;
  				 });
		      measures[s.simulation_num - 1] = on_arvs;
  		    })
  		  .allAgentsCreate(createDeterministicAgents)
  		  .timeAdjust("ARV_RATE", 1.0, 0, 1, sim::PROBABILITY)
  		  .parameter("ARV_RATE", {0.0, 0.0, 0.0, 0.0, 0.1} )
  		  .parameter("TIME_STEP", {0.5 } )
  		  .parameter("HIV_PREVALENCE", {0.0, 0.0, 0.0, 0.0} )
  		  .parameter("HIV_PREVALENCE_ARVS", {0.5} )
  		  .parameter("HIGH_RISK_PROPORTION", {0.5 } )
  		  .parameter("NUM_SIMULATIONS", {simulations})
		  .parameter("NUM_PARTNERSHIPS_LOW", {0} )
		  .parameter("NUM_PARTNERSHIPS_HIGH", {6} )).simulate();

  TESTRANGE(t, sim::mean(measures), 55, 65, "Mean of stochastic run.");

  std::cout << "Testing E8 year - WHO stage only to ARVs" << std::endl;
  measures.clear();
  measures.resize(simulations);
  sim::Simulation(sim::Options()
  		  .events({sim::advanceTimeEvent
  			, calcVariablesEvent
  			, arvsStochastic})
		  .beforeAllSimulations(setDefaultParameters)
  		  .afterEachSimulation([&measures](sim::Simulation &s) {
  		      unsigned on_arvs =
  			count_if(s.agents.begin(), s.agents.end(),
  				 [](const sim::Agent *a) {
  				   const HIVAgent *b = (HIVAgent *) a;
  				   return b->hiv == 5;
  				 });
		      measures[s.simulation_num - 1] = on_arvs;
  		    })
  		  .allAgentsCreate(createDeterministicAgents)
  		  .timeAdjust("ARV_RATE", 1.0, 0, 1, sim::PROBABILITY)
  		  .parameter("ARV_RATE", {0.1, 0.2, 0.3, 0.5, 0.0} )
  		  .parameter("TIME_STEP", { sim::YEAR } )
  		  .parameter("HIV_PREVALENCE", {0.1, 0.2, 0.3, 0.1} )
  		  .parameter("HIV_PREVALENCE_ARVS", {0.1} )
  		  .parameter("HIGH_RISK_PROPORTION", {0.5 } )
  		  .parameter("NUM_SIMULATIONS", {simulations})
		  .parameter("NUM_PARTNERSHIPS_LOW", {0} )
		  .parameter("NUM_PARTNERSHIPS_HIGH", {6} )).simulate();

  TESTRANGE(t, sim::mean(measures), 770, 792, "Mean of stochastic run.");

  std::cout << "Testing E9 - one year" << std::endl;
  measures.clear();
  measures.resize(simulations);
  sim::Simulation(sim::Options()
		  .events({sim::advanceTimeEvent
			, calcVariablesEvent
			, increasePopulationStochastic})
		  .beforeAllSimulations(setDefaultParameters)
  		  .afterEachSimulation([&measures](sim::Simulation &s) {
		      measures[s.simulation_num - 1] = s.agents.size();
  		    })
		  .agentCreate(create_hiv_agent)
		  .timeAdjust("GROWTH", 1.0, 0, 1, sim::PROBABILITY)
		  .parameter("HIGH_RISK_PROPORTION", {0.5 } )
		  .parameter("TIME_STEP", {1.0} )
		  .parameter("GROWTH", {0.01} )
		  .parameter("RISK_SUSCEPTIBLE", {0.0} )
		  .parameter("RISK_SUSCEPTIBLE_STDEV", {0.0} )
		  .parameter("RISK_UNINFECTIOUS", {0.0} )
		  .parameter("RISK_UNINFECTIOUS_STDEV", {0.0} )
		  .parameter("PARTNER_FORMATION_RATE", {0.0} )
		  .parameter("PARTNER_FORMATION_RATE_STDEV", {0.0} )
		  .parameter("HIGH_RISK_PROPORTION", {0.5 } )
		  .parameter("NUM_SIMULATIONS", {simulations})
		  .parameter("NUM_PARTNERSHIPS_LOW", {0} )
		  .parameter("NUM_PARTNERSHIPS_HIGH", {6} )
		  .parameter("PREVALENCE", {0.0} )).simulate();
  TESTRANGE(t, sim::mean(measures), 1205, 1235, "Mean of stochastic run.");

  std::cout << "Testing E9 - half year" << std::endl;
  measures.clear();
  measures.resize(simulations);
  sim::Simulation(sim::Options()
		  .events({sim::advanceTimeEvent
			, calcVariablesEvent
			, increasePopulationStochastic})
		  .beforeAllSimulations(setDefaultParameters)
  		  .afterEachSimulation([&measures](sim::Simulation &s) {
		      measures[s.simulation_num - 1] = s.agents.size();
  		    })
		  .agentCreate(create_hiv_agent)
		  .timeAdjust("GROWTH", 1.0, 0, 1, sim::PROBABILITY)
		  .parameter("HIGH_RISK_PROPORTION", {0.5 } )
		  .parameter("TIME_STEP", {0.5} )
		  .parameter("GROWTH", {0.01} )
		  .parameter("RISK_SUSCEPTIBLE", {0.0} )
		  .parameter("RISK_SUSCEPTIBLE_STDEV", {0.0} )
		  .parameter("RISK_UNINFECTIOUS", {0.0} )
		  .parameter("RISK_UNINFECTIOUS_STDEV", {0.0} )
		  .parameter("PARTNER_FORMATION_RATE", {0.0} )
		  .parameter("PARTNER_FORMATION_RATE_STDEV", {0.0} )
		  .parameter("HIGH_RISK_PROPORTION", {0.5 } )
		  .parameter("NUM_SIMULATIONS", {simulations})
		  .parameter("NUM_PARTNERSHIPS_LOW", {0} )
		  .parameter("NUM_PARTNERSHIPS_HIGH", {6} )
		  .parameter("PREVALENCE", {0.0} )).simulate();
  TESTRANGE(t, sim::mean(measures), 1205, 1235, "Mean of stochastic run.");

  // Test E10
  std::cout << "Testing E10 time step year" << std::endl;
  measures.clear();
  measures.resize(simulations);
  sim::Simulation(sim::Options()
		  .events({sim::advanceTimeEvent
			, calcVariablesEvent
			, diseaseProgressionStochastic})
		  .beforeAllSimulations(setDefaultParameters)
		  .afterEachSimulation([&measures](sim::Simulation &s) {
		      unsigned stage_1 =
			count_if(s.agents.begin(), s.agents.end(),
				 [](const sim::Agent *a) {
				   const HIVAgent *b = (HIVAgent *) a;
				   return b->hiv == 1;
				 });
		      measures[s.simulation_num - 1] = stage_1;
		    })
		  .agentCreate(create_hiv_agent)
		  .timeAdjust("DISEASE_CHANGE_RATE", 1.0, 0, 1,
			      sim::PROBABILITY)
		  .parameter("HIGH_RISK_PROPORTION", {0.5 } )
		  .parameter("DISEASE_CHANGE_RATE", {0.1, 0.2, 0.3} )
		  .parameter("TIME_STEP", {1.0 } )
		  .parameter("HIV_PREVALENCE", {1.0, 0.0, 0.0, 0.0} )
		  .parameter("NUM_YEARS", {4})
		  .parameter("NUM_SIMULATIONS", {simulations})
		  .parameter("NUM_PARTNERSHIPS_LOW", {0} )
		  .parameter("NUM_PARTNERSHIPS_HIGH", {6} )
		  .parameter("HIGH_RISK_PROPORTION", {0.5 } )).simulate();
  TESTRANGE(t, sim::mean(measures), 654, 659, "Mean of stochastic run.");

  std::cout << "Testing E10 time step half-year" << std::endl;
  measures.clear();
  measures.resize(simulations);
  sim::Simulation(sim::Options()
		  .events({sim::advanceTimeEvent
			, calcVariablesEvent
			, diseaseProgressionStochastic})
		  .beforeAllSimulations(setDefaultParameters)
		  .afterEachSimulation([&measures](sim::Simulation &s) {
		      unsigned stage_1 =
			count_if(s.agents.begin(), s.agents.end(),
				 [](const sim::Agent *a) {
				   const HIVAgent *b = (HIVAgent *) a;
				   return b->hiv == 1;
				 });
		      measures[s.simulation_num - 1] = stage_1;
		    })
		  .agentCreate(create_hiv_agent)
		  .timeAdjust("DISEASE_CHANGE_RATE", 1.0, 0, 1,
			      sim::PROBABILITY)
		  .parameter("HIGH_RISK_PROPORTION", {0.5 } )
		  .parameter("DISEASE_CHANGE_RATE", {0.1, 0.2, 0.3} )
		  .parameter("TIME_STEP", {0.5} )
		  .parameter("HIV_PREVALENCE", {1.0, 0.0, 0.0, 0.0} )
		  .parameter("NUM_YEARS", {4})
		  .parameter("NUM_SIMULATIONS", {simulations})
		  .parameter("NUM_PARTNERSHIPS_LOW", {0} )
		  .parameter("NUM_PARTNERSHIPS_HIGH", {6} )
		  .parameter("HIGH_RISK_PROPORTION", {0.5 } )).simulate();
  TESTRANGE(t, sim::mean(measures), 654, 659, "Mean of stochastic run.");

  std::cout << "Testing E12 time step one year, defaults" << std::endl;

  size_t fixed_agents;
  sim::Simulation(sim::Options()
		  .events({sim::advanceTimeEvent
			, calcVariablesEvent
			, InfectionRiskRandomMatching()})
		  .beforeAllSimulations(setDefaultParameters)
		  .afterEachSimulation([&](sim::Simulation &s) {
		      unsigned stage_1 =
			count_if(s.agents.begin(), s.agents.end(),
				 [](const sim::Agent *a) {
				   const HIVAgent *b = (HIVAgent *) a;
				   return b->hiv == 1;
				 });
		      fixed_agents = stage_1;
		    })
		  .agentCreate(create_hiv_agent)
		  .timeAdjust("NUM_PARTNERSHIPS_LOW", 1.0, 0, 1,
			      sim::LINEAR)
		  .timeAdjust("NUM_PARTNERSHIPS_HIGH", 1.0, 0, 1,
			      sim::LINEAR)
		  .parameter("HIGH_RISK_PROPORTION", {0.5 } )
		  .parameter("TIME_STEP", {1.0} )
		  .parameter("HIGH_RISK_PROPORTION", {0.5 } )
		  .parameter("INFECTION_PARAMETERS_TIME", {1.0})
		  .parameter("RISK_TRANSMISSION_ACT", {0.01} )
		  .parameter("NUM_PARTNERSHIPS_LOW", {0} )
		  .parameter("NUM_PARTNERSHIPS_HIGH", {6} )
		  .parameter("HIV_PREVALENCE", {0.1, 0.00, 0.00, 0.00} )
		  .parameter("NUM_SIMULATIONS", {1})
		  .parameter("NUM_PARTNERSHIPS_LOW", {0} )
		  .parameter("NUM_PARTNERSHIPS_HIGH", {6} )
		  .parameter("HIV_PREVALENCE_ARVS", {0.00} )).simulate();


  size_t poisson_agents;
  std::cout << "Testing E12 time step one year, poisson" << std::endl;
  sim::Simulation(sim::Options()
		  .events({sim::advanceTimeEvent
			, calcVariablesEvent
			, InfectionRiskRandomMatching(poisson_encounters)})
		  .beforeAllSimulations(setDefaultParameters)
		  .afterEachSimulation([&](sim::Simulation &s) {
		      unsigned stage_1 =
			count_if(s.agents.begin(), s.agents.end(),
				 [](const sim::Agent *a) {
				   const HIVAgent *b = (HIVAgent *) a;
				   return b->hiv == 1;
				 });
		      poisson_agents = stage_1;
		    })
		  .agentCreate(create_hiv_agent)
		  .timeAdjust("NUM_PARTNERSHIPS_LOW", 1.0, 0, 1,
			      sim::LINEAR)
		  .timeAdjust("NUM_PARTNERSHIPS_HIGH", 1.0, 0, 1,
			      sim::LINEAR)
		  .parameter("HIGH_RISK_PROPORTION", {0.5 } )
		  .parameter("TIME_STEP", {1.0} )
		  .parameter("HIGH_RISK_PROPORTION", {0.5 } )
		  .parameter("INFECTION_PARAMETERS_TIME", {1.0})
		  .parameter("RISK_TRANSMISSION_ACT", {0.01} )
		  .parameter("NUM_PARTNERSHIPS_LOW", {0} )
		  .parameter("NUM_PARTNERSHIPS_HIGH", {6} )
		  .parameter("HIV_PREVALENCE", {0.1, 0.00, 0.00, 0.00} )
		  .parameter("NUM_SIMULATIONS", {1})
		  .parameter("NUM_PARTNERSHIPS_LOW", {0} )
		  .parameter("NUM_PARTNERSHIPS_HIGH", {6} )
		  .parameter("HIV_PREVALENCE_ARVS", {0.00} )).simulate();
  TESTEQ(t, (fixed_agents < poisson_agents), true,
	 "More infections by poisson than fixed number of partners");


  size_t sticky_agents;
  std::cout << "Testing E12 time step one year, poisson, sticky" << std::endl;
  sim::Simulation(sim::Options()
		  .events({sim::advanceTimeEvent
			, calcVariablesEvent
			, InfectionRiskRandomMatching(poisson_encounters,
						      random_sticky_partner)})
		  .beforeAllSimulations(setDefaultParameters)
		  .afterEachSimulation([&](sim::Simulation &s) {
		      unsigned stage_1 =
			count_if(s.agents.begin(), s.agents.end(),
				 [](const sim::Agent *a) {
				   const HIVAgent *b = (HIVAgent *) a;
				   return b->hiv == 1;
				 });
		      sticky_agents = stage_1;
		    })
		  .agentCreate(create_hiv_agent)
		  .timeAdjust("NUM_PARTNERSHIPS_LOW", 1.0, 0, 1,
			      sim::LINEAR)
		  .timeAdjust("NUM_PARTNERSHIPS_HIGH", 1.0, 0, 1,
			      sim::LINEAR)
		  .parameter("HIGH_RISK_PROPORTION", {0.5 } )
		  .parameter("TIME_STEP", {1.0} )
		  .parameter("HIGH_RISK_PROPORTION", {0.5 } )
		  .parameter("INFECTION_PARAMETERS_TIME", {1.0})
		  .parameter("RISK_TRANSMISSION_ACT", {0.01} )
		  .parameter("NUM_PARTNERSHIPS_LOW", {0} )
		  .parameter("NUM_PARTNERSHIPS_HIGH", {6} )
		  .parameter("HIV_PREVALENCE", {0.1, 0.00, 0.00, 0.00} )
		  .parameter("NUM_SIMULATIONS", {1})
		  .parameter("NUM_PARTNERSHIPS_LOW", {0} )
		  .parameter("NUM_PARTNERSHIPS_HIGH", {6} )
		  .parameter("HIV_PREVALENCE_ARVS", {0.00} )).simulate();
  TESTEQ(t, (sticky_agents < poisson_agents), true,
	 "More infections without partner stickiness");


  std::cout << "Testing E12 time step one year, poisson, random k" << std::endl;
  size_t random_k_agents;
  sim::Simulation(sim::Options()
		  .events({sim::advanceTimeEvent
			, calcVariablesEvent
			, InfectionRiskRandomMatching
			(poisson_encounters,
			 SampleKPartners(700,
					 distance_simple))})
		  .beforeAllSimulations(setDefaultParameters)
		  .afterEachSimulation([&](sim::Simulation &s) {
		      unsigned stage_1 =
			count_if(s.agents.begin(), s.agents.end(),
				 [](const sim::Agent *a) {
				   const HIVAgent *b = (HIVAgent *) a;
				   return b->hiv == 1;
				 });
		      random_k_agents = stage_1;
		    })
		  .agentCreate(create_hiv_agent)
		  .timeAdjust("NUM_PARTNERSHIPS_LOW", 1.0, 0, 1,
			      sim::LINEAR)
		  .timeAdjust("NUM_PARTNERSHIPS_HIGH", 1.0, 0, 1,
			      sim::LINEAR)
		  .parameter("HIGH_RISK_PROPORTION", {0.5 } )
		  .parameter("TIME_STEP", {1.0} )
		  .parameter("HIGH_RISK_PROPORTION", {0.5 } )
		  .parameter("INFECTION_PARAMETERS_TIME", {1.0})
		  .parameter("RISK_TRANSMISSION_ACT", {0.01} )
		  .parameter("NUM_PARTNERSHIPS_LOW", {0} )
		  .parameter("NUM_PARTNERSHIPS_HIGH", {6} )
		  .parameter("HIV_PREVALENCE", {0.1, 0.00, 0.00, 0.00} )
		  .parameter("NUM_SIMULATIONS", {1})
		  .parameter("NUM_PARTNERSHIPS_LOW", {0} )
		  .parameter("NUM_PARTNERSHIPS_HIGH", {6} )
		  .parameter("HIV_PREVALENCE_ARVS", {0.00} )).simulate();
  TESTEQ(t, (sticky_agents < poisson_agents), true,
	 "More infections without partner stickiness");
  TESTRANGE(t, fabs( (double) random_k_agents - (double) sticky_agents),
	    0, 15, "No large difference between random_k and poisson");
  TESTRANGE(t, random_k_agents, 140, 170, "reasonable range for random_k");

  std::cout << "Testing E12 time step one week, poisson, random k" << std::endl;
  size_t random_k_week;
  sim::Simulation(sim::Options()
		  .events({sim::advanceTimeEvent
			, calcVariablesEvent
			, InfectionRiskRandomMatching
			(poisson_encounters,
			 SampleKPartners(700,
					 distance_simple))})
		  .beforeAllSimulations(setDefaultParameters)
		  .afterEachSimulation([&](sim::Simulation &s) {
		      unsigned stage_1 =
			count_if(s.agents.begin(), s.agents.end(),
				 [](const sim::Agent *a) {
				   const HIVAgent *b = (HIVAgent *) a;
				   return b->hiv == 1;
				 });
		      random_k_week = stage_1;
		    })
		  .agentCreate(create_hiv_agent)
		  .timeAdjust("NUM_PARTNERSHIPS_LOW", 1.0, 0, 1,
			      sim::LINEAR)
		  .timeAdjust("NUM_PARTNERSHIPS_HIGH", 1.0, 0, 1,
			      sim::LINEAR)
		  .parameter("HIGH_RISK_PROPORTION", {0.5 } )
		  .parameter("TIME_STEP", {sim::WEEK} )
		  .parameter("HIGH_RISK_PROPORTION", {0.5 } )
		  .parameter("INFECTION_PARAMETERS_TIME", {1.0})
		  .parameter("RISK_TRANSMISSION_ACT", {0.01} )
		  .parameter("NUM_PARTNERSHIPS_LOW", {0} )
		  .parameter("NUM_PARTNERSHIPS_HIGH", {6} )
		  .parameter("HIV_PREVALENCE", {0.1, 0.00, 0.00, 0.00} )
		  .parameter("NUM_SIMULATIONS", {1})
		  .parameter("NUM_PARTNERSHIPS_LOW", {0} )
		  .parameter("NUM_PARTNERSHIPS_HIGH", {6} )
		  .parameter("HIV_PREVALENCE_ARVS", {0.00} )).simulate();
  TESTRANGE(t, fabs( (double) random_k_agents - (double) random_k_week),
	    0, 10, "No large difference between week and year");

  std::cout << "Testing E12 time step one week, poisson, weighted" << std::endl;
  size_t weighted_agents;
  sim::Simulation(sim::Options()
		  .events({sim::advanceTimeEvent
			, calcVariablesEvent
			, InfectionRiskRandomMatching
			(poisson_encounters,
			 SampleKPartners(700,
					 distance_simple),
			 weighted_shuffle)})
		  .beforeAllSimulations(setDefaultParameters)
		  .afterEachSimulation([&](sim::Simulation &s) {
		      unsigned stage_1 =
			count_if(s.agents.begin(), s.agents.end(),
				 [](const sim::Agent *a) {
				   const HIVAgent *b = (HIVAgent *) a;
				   return b->hiv == 1;
				 });
		      weighted_agents = stage_1;
		    })
		  .agentCreate(create_hiv_agent)
		  .timeAdjust("NUM_PARTNERSHIPS_LOW", 1.0, 0, 1,
			      sim::LINEAR)
		  .timeAdjust("NUM_PARTNERSHIPS_HIGH", 1.0, 0, 1,
			      sim::LINEAR)
		  .parameter("HIGH_RISK_PROPORTION", {0.5 } )
		  .parameter("TIME_STEP", {sim::WEEK} )
		  .parameter("HIGH_RISK_PROPORTION", {0.5 } )
		  .parameter("INFECTION_PARAMETERS_TIME", {1.0})
		  .parameter("RISK_TRANSMISSION_ACT", {0.01} )
		  .parameter("NUM_PARTNERSHIPS_LOW", {0} )
		  .parameter("NUM_PARTNERSHIPS_HIGH", {6} )
		  .parameter("HIV_PREVALENCE", {0.1, 0.00, 0.00, 0.00} )
		  .parameter("NUM_SIMULATIONS", {1})
		  .parameter("NUM_PARTNERSHIPS_LOW", {0} )
		  .parameter("NUM_PARTNERSHIPS_HIGH", {6} )
		  .parameter("HIV_PREVALENCE_ARVS", {0.00} )).simulate();
  TESTRANGE(t, fabs( (double) random_k_agents - (double) weighted_agents),
	    0, 10, "No large deviation by weighted agents");

  std::cout << "Testing E12 time step one week, poisson, cluster" << std::endl;
  size_t cluster_agents;
  sim::Simulation(sim::Options()
		  .events({sim::advanceTimeEvent
			, calcVariablesEvent
			, InfectionRiskRandomMatching
			(poisson_encounters,
			 SampleKPartners(700,
					 distance_simple),
			 cluster_shuffle)})
		  .beforeAllSimulations(setDefaultParameters)
		  .afterEachSimulation([&](sim::Simulation &s) {
		      unsigned stage_1 =
			count_if(s.agents.begin(), s.agents.end(),
				 [](const sim::Agent *a) {
				   const HIVAgent *b = (HIVAgent *) a;
				   return b->hiv == 1;
				 });
		      cluster_agents = stage_1;
		    })
		  .agentCreate(create_hiv_agent)
		  .timeAdjust("NUM_PARTNERSHIPS_LOW", 1.0, 0, 1,
			      sim::LINEAR)
		  .timeAdjust("NUM_PARTNERSHIPS_HIGH", 1.0, 0, 1,
			      sim::LINEAR)
		  .parameter("HIGH_RISK_PROPORTION", {0.5 } )
		  .parameter("TIME_STEP", {sim::WEEK} )
		  .parameter("HIGH_RISK_PROPORTION", {0.5 } )
		  .parameter("INFECTION_PARAMETERS_TIME", {1.0})
		  .parameter("RISK_TRANSMISSION_ACT", {0.01} )
		  .parameter("NUM_PARTNERSHIPS_LOW", {0} )
		  .parameter("NUM_PARTNERSHIPS_HIGH", {6} )
		  .parameter("HIV_PREVALENCE", {0.1, 0.00, 0.00, 0.00} )
		  .parameter("NUM_SIMULATIONS", {1})
		  .parameter("CLUSTER_SIZE", {100})
		  .parameter("NUM_PARTNERSHIPS_LOW", {0} )
		  .parameter("NUM_PARTNERSHIPS_HIGH", {6} )
		  .parameter("HIV_PREVALENCE_ARVS", {0.00} )).simulate();
  TESTRANGE(t, fabs( (double) random_k_agents - (double) cluster_agents),
	    0, 10, "No large deviation by cluster agents");


  std::cout << "Testing E12 time step one week, poisson, reverse cluster"
	    << std::endl;
  size_t reverse_cluster_agents;
  sim::Simulation(sim::Options()
		  .events({sim::advanceTimeEvent
			, calcVariablesEvent
			, InfectionRiskRandomMatching
			(poisson_encounters,
			 SampleKPartners(700,
					 distance_simple),
			 cluster_shuffle, true)})
		  .beforeAllSimulations(setDefaultParameters)
		  .afterEachSimulation([&](sim::Simulation &s) {
		      unsigned stage_1 =
			count_if(s.agents.begin(), s.agents.end(),
				 [](const sim::Agent *a) {
				   const HIVAgent *b = (HIVAgent *) a;
				   return b->hiv == 1;
				 });
		      reverse_cluster_agents = stage_1;
		    })
		  .agentCreate(create_hiv_agent)
		  .timeAdjust("NUM_PARTNERSHIPS_LOW", 1.0, 0, 1,
			      sim::LINEAR)
		  .timeAdjust("NUM_PARTNERSHIPS_HIGH", 1.0, 0, 1,
			      sim::LINEAR)
		  .parameter("HIGH_RISK_PROPORTION", {0.5 } )
		  .parameter("TIME_STEP", {sim::WEEK} )
		  .parameter("HIGH_RISK_PROPORTION", {0.5 } )
		  .parameter("INFECTION_PARAMETERS_TIME", {1.0})
		  .parameter("RISK_TRANSMISSION_ACT", {0.01} )
		  .parameter("NUM_PARTNERSHIPS_LOW", {0} )
		  .parameter("NUM_PARTNERSHIPS_HIGH", {6} )
		  .parameter("HIV_PREVALENCE", {0.1, 0.00, 0.00, 0.00} )
		  .parameter("NUM_SIMULATIONS", {1})
		  .parameter("CLUSTER_SIZE", {100})
		  .parameter("NUM_PARTNERSHIPS_LOW", {0} )
		  .parameter("NUM_PARTNERSHIPS_HIGH", {6} )
		  .parameter("HIV_PREVALENCE_ARVS", {0.00} )).simulate();
  TESTEQ(t, (reverse_cluster_agents < cluster_agents), true,
	 "Cluster agents sorted by least risky results in fewer infections");

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
