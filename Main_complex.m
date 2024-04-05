close all; clear all; clc
%This is the main-script for a complex valued watertank EIT
%simulation without inverse crime.

%% Create simulated data
disp('Creating simulated data');

simuname = 'complexsimu';
nel = 32;
vincl = true(nel,nel);
vincl = vincl(:);

load_data = 0;%Do we run the simulation or just load the data.
if load_data
    load([simuname 'data.mat']);
    disp(['Loaded simulated data from file ' simuname 'data.mat']);
else
    load('meshfiles/KIT_mesh_dense.mat');%Load mesh variables
    ginvsimu = ginv(:,1:2);
    Hinvsimu = Hinv;
    g1stsimu = g;
    H1stsimu = H;
    elfaces1stsimu = elfaces1st;
    sigmasimu1 = GenerateEllipseInGrid(ginvsimu, Hinvsimu, 1, 10, 4e-2, 3e-2, 5e-2, 0e-2, 1e-2, 1);%Generate a single blob as a target
    sigmasimu2 = GenerateEllipseInGrid(ginvsimu, Hinvsimu, 1, 10, 4e-2, 3e-2, -5e-2, 0e-2, 1e-2, 1);%Generate a single blob as a target
    sigmasimu = Estimate_vec({sigmasimu1; sigmasimu2});
    z = 1e-9*ones(nel,1);
    fmsimu = ForwardMesh1st(g1stsimu, H1stsimu, elfaces1stsimu);
    imesh.g = ginvsimu;
    imesh.H = Hinvsimu;
    fmsimu.SetInverseMesh(imesh);
    simusolver = EITFEM_complex(fmsimu);
    simusolver.zeta = z;
    [Umeas, Imeas] = Simulation(fmsimu, [], sigmasimu, [], simusolver, [], [0 0 0 0]);
    save([simuname 'data.mat'], 'Umeas', 'Imeas', 'sigmasimu', 'z', 'ginvsimu', 'Hinvsimu');
    disp(['simulation run succesfully, results saved in ' simuname 'data.mat']);
end

%% Load meshes for inverse problem
%Load meshes for solving the inverse problem, i.e. a 3D forward mesh and a
%2D inverse mesh.
load('meshfiles/KIT_mesh.mat');
ginv = ginv(:,1:2);
g1st = g;
H1st = H;
fm = ForwardMesh1st(g1st, H1st, elfaces1st);
imesh.g = ginv;
imesh.H = Hinv;
fm.SetInverseMesh(imesh);
disp('meshfiles loaded');


%% Setting up the inverse solver

%Set up the forward problem solver:
solver = EITFEM_complex(fm);
solver.sigmamin = 1e-9;
solver.zeta = 1e-9*ones(nel,1);
solver.mode = 'potential';
solver.vincl = vincl;
solver.Iel = Imeas;
solver.Uel = Umeas;
solver.SetInvGamma(1e-4, 3e-2);
disp('Forward problem solver set up')

%Set up the Total Variation prior:
ec = 10;%expected change
TVPrior = PriorTotalVariation(ginv, Hinv, [ec; ec]);
TVPrior.sigmaind = [1; 2];%apply total variation to estimates with indices 1 and 2
disp('TV prior set up');

%Set up the positivity prior:
PosiPrior = PriorPositivityParabolic([1e-5; 1e-5], [1e2; 1e2]);
disp('Positivity prior set up');

%Set up the plotter:
plotter = Plotter(ginv, Hinv);
plotter.plottypes = [1; 1];

%Finally, set up the inverse problem solver, in this case GN with
%linesearch:
resobj = cell(3, 1);
resobj{1} = solver;
resobj{2} = TVPrior;
resobj{3} = PosiPrior;
InvSolver = SolverLinesearch(resobj);
InvSolver.maxIter = 100;
InvSolver.Plotter = plotter;

%Make the initial guess and start solving!
sigmainitial = ones(size(ginv,1),1);
esti = Estimate_vec({sigmainitial; sigmainitial});
disp('All set! Beginning to solve the inverse problem.')
reco = InvSolver.solve(esti);


