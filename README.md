# Ehrenfest-like particle system in C++

## Brief summary
This repository contains a setup for a ballistic particle model in two urns. It extends the Ehrenfest model and has been used in [this][1] and [this][2] publication.

The goal of this code is to explore long-term behaviour of the model with different urn transition mechanisms. Two have been implemented, and the code is flexible enough to support other custom implementations.
Rather than a time-discretized ODE system to simulate particle trajectories, it uses a discrete-event simulation that predicts collisions and interpolates particle positions between these collisions.

The code has some scripts that perform pre-processing jobs (creating parameter files to explore parameter regimes, setting up batch jobs that can run in parallel) and post-processing jobs (visualizing trajectories, creating graphs and images used in the publications).

## Requirements
The source code for the simulation is written in C++ (using some c++14 features) and has no external dependencies. The test suite relies on [Boost][3].
The code can be compiled with [CMake][4]. If not installed, get it via your favorite package manager.
Installation and executation has been tested on Linux and Mac, and perhaps also works on Windows.

## Installation instructions
Installation is pretty straightforward: clone the repository and run CMake and make to obtain the executables. The following workflow suffices:
```bash
git clone https://github.com/0mar/particular.git
cd particular
## build the executables
mkdir build
cd build
cmake -DCMAKE_BUILD_TYPE=Release ..
make
```
For development and debugging (getting helpful error messages when stuff goes wrong, but a significant slower execution), replace `Release` with `Debug`.

## Running the code
This framework has 5 executables:
 - `particular`, which mainly functions as a demonstration of the model
 - `examinations`, featuring detailed data collections for the anti-diffusion phenomenon
 - `single_channel`, containing functions to run the experiments and collect the data in [this][1] paper
 - `double_channel`, containing functions to run the experiments and collect the data in [this][2] paper
 - `test_particular`, to run the unit test suite

The last two executables take a list of parameters as arguments. These are not really meant to be run manually.
The Python scripts `create_single_channel_batch.py` and `create_double_channel_batch.py` respectively create parameter files that these executables accept.

Furthermore, there are various scripts that post-process the data from these executables:
 - `visualise.py` creates a TKinter animation that one can use to see the phenomenon in action
 - `plot_data.py` creates the plots used in the papers
 - `plot_thermalisation.py` creates figures that illustrate long-term behavior

## Unit test suite
This code has a test suite that mainly tests the computational geometry and the custom data structures. It relies on Boost and is run after `make test` with `./test_particular`. If you make any changes to the geometry, adding corresponding tests is highly recommended. Computational arithmetic is fickle and round-off errors quickly accumulate, especially when dealing with many particles and long run-times.

[1]: https://doi.org/10.1088/1751-8121/ab94ec
[2]: https://doi.org/10.1103/PhysRevE.103.032119
[3]: https://www.boost.org/
[4]: https://cmake.org/
