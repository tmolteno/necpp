#RANT: Ruby Wire Antenna Optimiser

Author:: Tim Molteno (tim@physics.otago.ac.nz)
Copyright:: Copyright (C) Tim Molteno 2008-2013.
License:: Licensed under the GNU GPL v3

##Overview

This software uses genetic programming to optimize wire antenna structures. The genetic
programming component of the software is written in the Ruby programming language, and
the simulations are done using the Ruby binding to the nec2++ antenna simulation
software.

This software contains two components a Gp server object, and a SimulationClient component that does
antenna fitness evaluations.

The main GP algorithm is described by the Gp object
that places Individual members of the Population into an EvalQueue. The EvalQueue is serviced
by 


##Howto

To modify for a different fitness, The methods that need to be overridden are:
* Antenna.get_statistics(freq)
* Individual.fitness_calc 
* Individual.fitness_from_stats(freq)

First issue the command 'make test10000'. Then on one or more machines, issue the command 'make client'. 
The computation client (there can be hundreds of them) are essentially grid compute engines that simulate
antenna using nec2++.


##ChangeLog

* 2013/03/19: Cleaned up code to re-segment wires.
