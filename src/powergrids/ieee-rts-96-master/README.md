# IEEE RTS-96

The IEEE Reliability Test System-1996 is a benchmark system for power grids. It consists of 73 buses and 108 lines with 33 generators.

This repository is supposed to provide code for building dynamical test cases with PowerDynamics.jl and NetworkDynamics.jl.

In the long run this project will eventually turn into a package that provides Matpower / PSSE source files of different benchmark systems and uses the PowerModels.jl parsers for building dynamical testcases.

![Graph structure of the IEEE RTS 1996](rts-96_graph.png)

## Scripts

1. *ND_Script.jl* is an example script that shows how to simulate the system with the NetworkDynamics.jl package. At the current state generator nodes are modeled as swing equations. Load nodes are modeled as algebraic contraints or dynamically, following Bergen & Hill (1981). In *ND_helper_functions.jl* you can find some convenience functions for parsing and plotting the testcase.
2. *PowerModels_Script.jl* is an example script that shows how to parse the system from a .raw file, calculates the AC power flow and exports the results to *Flow.csv*

**Note: This is work under progress. Some details have not been implemented yet!**

## Data

In *Data* I provide .csv files for conveniently building the testcases in Julia.

| File | Data |
|------|------|
| *Bus.csv* | Number, ID, Names, Type & Base Voltage of Buses |
| *Line.csv* | Source / Destination Buses, Reactances & Resistances of Lines |
| *Load.csv* | Bus ID, Active & Reactive Power of Loads |
| *Generator.csv* | Bus ID, Voltage Setpoints and Active Power of Generators |
| *Flow.csv* | Voltage Angles and Magnitudes of Power Flow Solution |

Other Data Formats:
* *Data/Original_Data* contains the tables of the original publication in a text format
* *Data/TAMU_Data* contains a .raw file of the testcase
* *Data/NICTA_Data* contains an .m file of the testcase

## Literature

In *Literature* you can find the original IEEE paper describing the RTS-1996 testcase.

## To Do

Here is a list of future improvements:

* detailed generator modeling (3rd/4th order, heterogeneous parameters)
* include shunts (pi-model lines)
* build Matpower testcase and parse it with PowerModels.jl
* include spatial data of the buses
* restructure repository as a Julia package that can directly be used in a script

I wrote GitLab issues for all of these future improvements.