/*
 * Used only for testing the Simulator.
 */

#include "sim.hh"

int main(int argc, char **argv)
{
  sim::Simulation(sim::Options()
		  .events({sim::advanceTimeEvent})
		  .commandLine(argc, argv)).simulate();
  return 0;
}
