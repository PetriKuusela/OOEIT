# OOEIT
OOEIT - Object Oriented Electrical Impedance Tomography, an open source code package for Matlab.
## Purpose
The OOEIT package is aimed for anyone who wants to compute EIT reconstructions without having to understand the functioning of the whole reconstruction algorithm, as well as for anyone with deeper understanding who wants to create their own reconstruction codes.
## Getting started
The code package is ready for use right after downloading. Before using OOEIT, each Matlab session should run the initialize_OOEIT(); function, which adds the relevant subfolders to the search path. To get an idea of the usage you can check the main files included, which provide a starting point for few of the possible basic usage scenarios. Some explanations can also be found in the "Usage" section of this README-file. For a description of the structure of the OOEIT code package see the "Structure" section.
## Usage
The OOEIT package can be used both for simulating EIT measurements and computing EIT reconstructions.
### Simulating EIT measurements
To simulate EIT measurements, the user has to initialize a forward problem solver (FPS). There are a few different FPSs included in OOEIT, but the basic one is the EITFEM class. To initialize an EITFEM object, a mesh object is required. The basic mesh classes in OOEIT to use are the ForwardMesh1st and ForwardMesh2nd classes. Their initializers require the mesh nodes, element connectivity matrices and a cell array containing elements facing each electrode. The FPS is ready for use right after initialization, and measurements on the given mesh can be simulated by SolveForwardVec(sigma) function, where sigma is a vector containing the conductivity values at the mesh nodes.
Although all the other necessary parameters than the mesh are set by the initializer, the user may want to modify several of them:
- mode: Set to 'current' for current injection device simulations and 'potential' for potential injection device.
- Mpat: Measurement pattern. If your potential measurements are not agains a common ground, this matrix defines the pattern of measuring potentials
- Iel or Uel: The electrode currents or potentials injected through the target, i.e. the injection pattern
- zeta: A vector containing the contact impedances
- sigmamin: The minimum value of conductivity (default 10<sup>-6</sup>). All node values below this will be replaced by this value, so if your conductivity consists of smaller values than this, you should lower sigmamin accordingly.
### Computing EIT Reconstructions
The EIT reconstruction is computed by an inverse solver. The basic inverse solver class in OOEIT is the SolverGN, which solves the presented optimization problem by Gauss-Newton method. Initializing a SolverLinesearch object requires a cell array of objective functionals, the sum of which the object tries to minimize. In other words, this cell array contains the FPS and the prior objects, which have to be initialized before initializing the inverse Solver. The data required for initializing the priors depends on which prior classes are used. However, their initialization should be very straight forward, and by copying behaviour from included main files, the user should be able to easily use the prior classes.
After initializing the inverse solver, the reconstruction can be solved with the solve(sigma_initial) function, where sigma_initial is the initial guess, i.e. the starting point of the optimization algorithm. There are some properties the user may wish to set before this, however. Here is a list of the most important ones:
- Plot flags: plotLinesearch, plotIterations, plotData, plotConvergence and plotUQ. These control which plots are to be plotted during the reconstruction process.
- estop and nStepsBack: These control the stopping criterion of the Gauss-Newton algorithm. More precisely, the algorithm stops at iteration i, whenever the iteration i - nStepsBack has the sum of optimization functions less than the sum of optimization functions at i plus estop.
- nstd: If uncertainty quantification (UQ) is plotted, this defines the number of standard deviations used as its width.
- maxIter: Maximum number of Gauss-Newton iterations.
- maxIterInLine: Maximum number if iterations while trying to find the minimum along search direction inside each Gauss-Newton iteration.
- Plotter: A plotter object which can plot the estimate or UQ between iterations
- showSplitVals: A flag determining whether to plot each optimization function value separately during line search.
- TrueEst: The true value of the estimate, so that in can be compared to in the UQ lineplot.
## Structure
This will be updated later.
## Features
Currently, in addition to the very basic EIT reconstruction, OOEIT includes the following main features:
- Contact impedance estimation
- Complex EIT
- Non-linear difference reconstruction
- Support for 1st and 2nd order meshes
The following priors are included:
- Total variation prior
- Smoothness prior
- Positivity constraint
## Acknowledgements
The codes are based on earlier versions of EIT codes circulating at the University of Eastern Finland. Furhter, this work has been supported by the Research Council of Finland (the Finnish Center of Excellence of Inverse Modeling and Imaging, project number 353087), and Research Council of Finland's Flagship of Advanced Mathematics for Sensing, Imaging and Modelling, grant number 358944.
