close all; clear all; clc
%This is the main-script for a watertank EIT
%simulation where also the contact impedances are estimated without inverse
%crime.

initialize_OOEIT();%Make sure all the relevant paths are added

%% Create simulated data
disp('Creating simulated data');

load('meshfiles/Mesh2D_dense.mat');%Load mesh variables
fmsimu = ForwardMesh1st(g, H, elfaces);

sigmasimu = GenerateEllipse(g, 1, 10, 4e-2, 3e-2, 5e-2, 0e-2, 1e-2);%Generate a single blob as a target
z = 1e-6*ones(length(elfaces),1);%Set the contact impedances
z(3) = 1e-3;%A couple of higher ones to find through estimation
z(13) = 1e-2;
[Umeas, Imeas] = Simulation(fmsimu, [], sigmasimu, z, [], [], [0 0 0 0]);

%% Load meshes for inverse problem
%Load meshes for solving the inverse problem, i.e. a 3D forward mesh and a
%2D inverse mesh.
load('meshfiles/Mesh2D_sparse.mat');
fm = ForwardMesh1st(g, H, elfaces);
disp('meshfiles loaded');


%% Setting up the inverse solver

%Set up the forward problem solver:
solver = EITFEM(fm);
solver.Iel = Imeas;
solver.Uel = Umeas;
solver.SetInvGamma(1e-3, 3e-2);
solver.zind = 2;%The index of contact impedance estimate in the Estimate_vec object
disp('Forward problem solver set up')

%Set up the Total Variation prior:
%Start by calculating the value for alpha based on Gerardos method
ec = 1;%expected change
TVPrior = PriorTotalVariation(g, H, ec); %By default, TV is applied to estimate with index 1, which is the conductivity now
disp('TV prior set up');

%Set up the positivity prior:
PosiPrior = PriorPositivityParabolic([1e-5; 1e-9], [100, 1]);%Two values for each argument, for conductivity and for z
disp('Positivity prior set up');

%Set up the plotter:
plotter = Plotter(g, H);%a normal plotter for conductivity
plotterz = PlotterSimple();%and a plotter for z. This just plots the data in a lineplot
ph = PlotterHolder({plotter; plotterz});%The two plotters are then combined to a single PlotterHolder object

%Finally, set up the inverse problem solver, in this case GN with
%linesearch:
resobj = cell(3, 1);
resobj{1} = solver;
resobj{2} = TVPrior;
resobj{3} = PosiPrior;
InvSolver = SolverGN(resobj);
InvSolver.Plotter = ph;

%Make the initial guess and start solving!
sigmainitial = ones(size(g,1),1);
zinitial = 1e-3*ones(length(elfaces),1);
est = Estimate_vec({sigmainitial; zinitial});%This is the full initial guess

disp('All set! Beginning to solve the inverse problem.')
reco = InvSolver.solve(est);


