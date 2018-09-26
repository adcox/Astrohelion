![Build Status](https://travis-ci.org/adcox/astrohelion.svg?branch=master)

## Welcome!

Astrohelion (Astro - heel - eon) is a set of C++ tools that facilitiate 
multi-body trajectory design. Specifically, analysis in the circular restricted 
three-body problem (CR3BP), the CR3BP with low-thrust (CR3BP-LT), and bi-circular 
restricted four-body problem (BC4BP) is supported.

## Capabilities

So, what can Astrohelion do? There are several key functionalities:


* **Numerical Integration** Astrohelion leverages the ordinary differential 
	equation (ODE) solvers in the GSL library to numerically integrate the 
	equations of motion associated with a few common multi-body astrodynamical 
	systems. This propagation is controlled by a `SimEngine` object that 
	possesses many of the same features as Matlab's ODE solvers, such as event 
	detection. The engine creates a trajectory object that can then be operated 
	on for further analysis.

* **Differential Corrections** Continuous trajectories may be discretized into a 
	series of nodes, and constraints can be applied to the nodes, segments 
	between nodes, and the trajectory as a whole. A multiple shooting 
	differential corrections algorithm is then applied to the constrained 
	"arcset" with the goal of satisfying all constraints. The algorithm itself, 
	contained in the `MultShootEngine` object, is problem-independent and can be
	applied to CR3BP, CR3BP-LT, or BC4BP arcsets using the same set of commands.

* **Continuation** Once a desirable trajectory has been computed (and possibly 
	corrected to meet certain constraints), continuation algorithms can be 
	applied to generate a continuous family of similar solutions. Both natural 
	parameter continuation (see `NatParamEngine`) and psuedo-arclength 
	continuation (see `PseudoArcEngine`) are avaiable, both of which are 
	agnostic to dynamical model.

* **Export to File** Trajectories are easily saved to file in MATLAB binary 
	format for further analysis or visualization.

## Limitations

This software has been developed to support MS and PhD research. Thus, while the
tools are fully operational and have been tested fairly well, they are not as
flexible or extensible as they could be. Please submit feedback via issue reports
if you encounter problems so that we can continue to improve this software! 
However, note that we make no guarantee that we'll be able to offer support; this
software is provided as is, without any quality guarnatee.

## Installation

The following steps describe a general installation process.

1. Copy the `installDeps.sh` script into the main directory from the 
	`build` directory.
1. Open installDeps.sh and set the `INSTALL_DIR`, `CONFIG_DIR`, and `SPICE_DIR` 
	variables to point to locations you have write access to.
1. Additionally, set the `SPK` and `TLS` variables to the names of SPICE kernels 
	you want to download to use for ephemeris data.
1. Run the installation script, i.e., `sh installDeps.sh`. This will download, 
	build, and install the libraries required to run Astrohelion
1. Copy `makefile` out of the `build` directory into the main directory and 
	adjust the `INSTALL_DIR` variable to point to the same directory as the 
	`INSTALL_DIR` variable in `installDeps.sh`
1. Compile the project source, i.e., run `make`
1. Install the Astrohelion libriaries and headers, i.e., run `make install`

## Documentation

Astrohelion documentation is available at 
[http://adcox.github.io/astrohelion/](http://adcox.github.io/astrohelion/). To 
build the latest documentation (which is stored in `docs/api/html`), run the command 
`make docs`. Document generation requires the `doxygen` package.
