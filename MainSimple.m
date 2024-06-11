%This main-file features the example codes presented in the OOEIT article
%(not yet published by the time of this comment). This is a very basic
%usecase of first simulating the EIT measurements and then reconstructing
%the conductivity distribution based on the measurements.
% 
% Note that in this example, the inverse crime is commited! That is, the
% same model (specifically, the same mesh as well) is used both for
% simulating the data and as the forward model while solving the inverse
% problem. For brief reasoning why it is bad, see the OOEIT article.
clear all;
close all;

InitializeOOEIT();%Make sure all the relevant paths are added to the search path

%Start by loading one of the predefined meshes
load('meshfiles/Mesh2D_dense.mat');
%And create a mesh-object from the loaded variables g, H and elfaces
mesh = ForwardMesh1st(g, H, elfaces);

%Next, initialize the forward problem solver (FPS). This is the object, that
%(utilizing the mesh-object given to it) computes forward problem solutions
%using the FEM approximation of the complete electrode model (CEM).
solver = EITFEM(mesh);

%The FPS options can be set by simple assignments. As an example, let us
%change the injection mode to potential injection. If you want to use
%current injection, it is used as default, and you can just comment this
%line out. Then you just have to make sure you assign the simulated
%measurements to the correct variable (solver.Uel instead of solver.Iel).
solver.mode = 'potential';

%Next, let us create the target. This auxiliary function creates a
%conductivity distribution featuring a single ellipse on a homogeneous
%background
sigma = GenerateEllipse(g, 1, 10, 3e-2, 4e-2, 3e-2, 0, 1e-2);

%Then we can compute the forward problem, which means simulating the
%measurements, by the following line (note that if
% solver.mode == 'current', the solution is the potentials, and if
% solver.mode == 'potential', the solution is the currents.)
Imeas = solver.SolveForwardVec(sigma);

%For the inverse problem solution, regularization is required, and it is
%implement in the form of priors. Let us next initialize a total variation
%(TV) prior. There is no hard-coded limit for the number of priors you can
%use simultaneously, and often it is desirable to have, for example, at
%least a positivity constraint (e.g. PriorPositivityParabolic) added as well.
TVPrior = PriorTotalVariation(g, H, 3);

%To initialize the inverse problem solver (or inverse solver, for short), 
%we need a cell array containing our optimized functionals, i.e. the FPS
%and the priors.
ofuns = {solver; TVPrior};
%Then initialize the inverse problem solver
invSolver = SolverGN(ofuns);

%The inverse solver has also many tunable parameters, which may be
%mdified by simply assigning new values; here we decrease the maximum
%number of iteration inside linesearch:
invSolver.maxIterInLine = 15;

%To plot the estimate, we initialize a plotter class, and assign it to the
%inverse solver
plotter = Plotter(g, H);
invSolver.plotter = plotter;

%You can plot the true sigma using the plotter class by calling:
%plotter.plot(sigma);

%Before starting to solve the inverse problem, make sure that the
%measurement data has been given to the FPS. Furthermore, if you change the
%solver mode to current injection, the measurement data has to be passed to
%solver.Uel, as the measurements will then be the electrode potentials.
solver.Iel = Imeas;

%Make an initial guess
sigmainitial = ones(size(g,1),1);
%And start solving the inverse problem.
reconstruction = invSolver.Solve(sigmainitial);
%The result will be a vector of
%conductivity values at the node points of the mesh, and can be plotted
%using the same plotter class given to the inverse solver by calling
%plotter.plot(reconstruction);

