# Piston

This little program simulates an ideal gas in a one dimensional container with a movable lid. The lid pushed down by gravitation only. It illustrates oscillations around thermal equilibrium. This only works, however, if the frequency of these oscillations is slow compared to the relaxation time of the system. Otherwise, the system does not stay in equilibrium and the fluctuations become much more complicated

The units are such that the mass of the lid is 1.

## How to run

The simulation is written in C because the problem doesn’t map to python well. Build the library with `make`.

Then you can run `python plot.py`.
