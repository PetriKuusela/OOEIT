close all; clear all; clc
%This is the main-script for a complex valued watertank EIT
%simulation without inverse crime.

initialize_OOEIT();%Make sure all the relevant paths are added

%% Create simulated data
disp('Creating simulated data');

load('meshfiles/Mesh2D_dense.mat');%Load mesh variables
fmsimu = ForwardMesh1st(g, H, elfaces);

sigmasimu1 = GenerateEllipse(g, 1, 10, 4e-2, 3e-2, 5e-2, 0e-2, 1e-2);%Generate a single blob as a target conductivity
sigmasimu2 = GenerateEllipse(g, 1, 10, 4e-2, 3e-2, -5e-2, 0e-2, 1e-2);%Generate a single blob as a target permittivity
sigmasimu = Estimate_vec({sigmasimu1; sigmasimu2});

%Create the forward problem solver which solves the complex FEM
simusolver = EITFEM_complex(fmsimu);
%Run the simulations measuring also the phase data of the signals
[Umeas, Imeas] = Simulation(fmsimu, [], sigmasimu, [], simusolver, [], [0 0 0 0]);

%% Load meshes for inverse problem
%Load meshes for solving the inverse problem, i.e. a 3D forward mesh and a
%2D inverse mesh.
load('meshfiles/Mesh2D_sparse.mat');
fm = ForwardMesh1st(g, H, elfaces);
disp('meshfiles loaded');


%% Setting up the inverse solver

%Set up the forward problem solver:
solver = EITFEM_complex(fm);
solver.Iel = Imeas;
solver.Uel = Umeas;
solver.SetInvGamma(1e-4, 3e-2);%Set the noise level (this affects how the difference of the forward model and measurements are weighed while optimizing)
disp('Forward problem solver set up')

%Set up the Total Variation prior:
ec = 3;%expected change
TVPrior = PriorTotalVariation(g, H, [ec; ec]);%here "ec" is duplicated, as the prior applies to both initial and final state conductivity
TVPrior.sigmaind = [1; 2];%apply total variation to estimates with indices 1 and 2
disp('TV prior set up');

%Set up the positivity prior:
PosiPrior = PriorPositivityParabolic([1e-5; 1e-5], [1e2; 1e2]);%two parameters for both arguments are given; one for each estimate
disp('Positivity prior set up');

%Set up the plotters:
plotter = Plotter(g, H);%Both estimates require their own plotter
plotter2 = Plotter(g, H);
plotter2.title = 'Permittivity';
ph = PlotterHolder({plotter; plotter2});%The two plotters are then combined to a single PlotterHolder object

%Finally, set up the inverse problem solver, in this case GN with
%linesearch:
resobj = cell(3, 1);%Collect all the functionals to be minimized into a single cell array
resobj{1} = solver;
resobj{2} = TVPrior;
resobj{3} = PosiPrior;
InvSolver = SolverGN(resobj);%Create the inverse solver object
InvSolver.Plotter = ph;%Set the plotter of the inverse solver

%Make the initial guess and start solving!
sigmainitial = ones(size(g,1),1);
esti = Estimate_vec({sigmainitial; sigmainitial});%The initial guess for both conductivity and permittivity is just ones
disp('All set! Beginning to solve the inverse problem.')
reco = InvSolver.solve(esti);


