close all; clear all; clc
%This is the driver-script for a watertank EIT
%simulation without inverse crime, computing the reconstruction using the
%non-linear difference method.

%% Create simulated data
disp('Creating simulated data');

load('meshfiles/Mesh2D_dense.mat');%Load mesh variables
fmsimu = ForwardMesh1st(g, H, elfaces);%And create a mesh object

sigmasimu_i = GenerateEllipse(g, 1, 2, 4e-2, 3e-2, 5e-2, 0e-2, 1e-2);%Generate a single blob as initial conductivity
sigmasimu = GenerateEllipse(g, 0, 10, 4e-2, 3e-2, -5e-2, 0e-2, 1e-2);%Generate another blob for the change in condcutivity
sigmasimu = sigmasimu + sigmasimu_i;%This is the final conductivity
z = 1e-6*ones(length(elfaces),1);%contact impedances
noise = [0 0 0 0];%A four-component noise model is used; see Simulation.m for details. Here we have no noise

%Then simulate the measurements:
[Umeas, Imeas, Umeas_i, Imeas_i] = Simulation(fmsimu, sigmasimu_i, sigmasimu, z, [], [], noise);

%% Load meshes for inverse problem
%Load meshes for solving the inverse problem. Here no separate mesh for the
%inverse problem is used.
load('meshfiles/Mesh2D_sparse.mat');
fm = ForwardMesh1st(g, H, elfaces);%Create a mesh object from the loaded variables
disp('meshfiles loaded');


%% Setting up the inverse solver

%Set up the forward problem solver, which uses the initial state measurements
solver = EITFEM(fm);
solver.Iel = Imeas_i;
solver.Uel = Umeas_i;
solver.SetInvGamma(1e-3, 3e-2);%Set the noise level (this affects how the difference of the forward model and measurements are weighed while optimizing)

%Set up the forward problem solver, which uses the final state measurements
solver_nld = EITFEM_NLD(fm);
solver_nld.Iel = Imeas;
solver_nld.Uel = Umeas;
solver_nld.SetInvGamma(1e-3, 3e-2);
disp('Forward problem solver set up')

%Set up the Total Variation prior:
ec = 10;%expected change
TVPrior = PriorTotalVariation(g, H, [ec; ec]);%here "ec" is duplicated, as the prior applies to both initial and final state conductivity
TVPrior.sigmaind = [1; 2];%These are the indices (in the Estimate_vec) of the estimates this prior applies to (here both of the estimates)
disp('TV prior set up');

%Set up the positivity prior:
PosiPrior = PriorPositivityParabolic(1e-5, 100);
PosiPrior.omitind = 2;%Tell the positivity prior that estimate number 2 is not to be penalized, as the change may be negative as well
disp('Positivity prior set up');

%Set up the plotter:
plotter = Plotter(g, H);%Both estimates require their own plotter
plotter.title = 'Initial conductivity';
plotter2 = Plotter(g, H);
plotter2.title = 'Conductivity change';
ph = PlotterHolder({plotter; plotter2});%The two plotters are then combined to a single PlotterHolder object

%Finally, set up the inverse problem solver, in this case GN with
%linesearch:
resobj = cell(4, 1);%Collect all the functionals to be minimized into a single cell array
resobj{1} = solver;
resobj{2} = solver_nld;
resobj{3} = TVPrior;
resobj{4} = PosiPrior;
InvSolver = SolverGN(resobj);%Create the inverse solver object
InvSolver.Plotter = ph;%Set the plotter of the inverse solver

%Make the initial guess and start solving!
sigmainitial = ones(size(g,1),1);%initial guess for the initial state of the conductivity
changeinitial = zeros(size(g,1),1);%initial guess for the change of conductivity
%Use the Estimate_vec class when having multiple parameters (here initial
%conductivity and the change) to estimate:
sigmainitial = Estimate_vec({sigmainitial; 0*sigmainitial});

disp('All set! Beginning to solve the inverse problem.')
reco = InvSolver.solve(sigmainitial);%Start solving the inverse problem


