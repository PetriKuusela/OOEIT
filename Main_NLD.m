close all; clear all; clc
%This is the driver-script for awatertank EIT
%simulation without inverse crime, computing the reconstruction using the
%non-linear difference method.

%% Create simulated data
disp('Creating simulated data');

simuname = 'NLDsimu';
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
    sigmasimu_i = GenerateEllipseInGrid(ginvsimu, Hinvsimu, 1, 2, 4e-2, 3e-2, 5e-2, 0e-2, 1e-2, 1);%Generate a single blob as a target
    sigmasimu = GenerateEllipseInGrid(ginvsimu, Hinvsimu, 0, 10, 4e-2, 3e-2, -5e-2, 0e-2, 1e-2, 1);
    sigmasimu = sigmasimu + sigmasimu_i;
    z = 1e-9*ones(nel,1);
    fmsimu = ForwardMesh1st(g1stsimu, H1stsimu, elfaces1stsimu);
    imesh.g = ginvsimu;
    imesh.H = Hinvsimu;
    fmsimu.SetInverseMesh(imesh);
    [Umeas, Imeas, Umeas_i, Imeas_i] = Simulation(fmsimu, sigmasimu_i, sigmasimu, z, [], [], [0 0 0 0]);
    save([simuname 'data.mat'], 'Umeas', 'Imeas', 'Umeas_i', 'Imeas_i', 'sigmasimu', 'z', 'ginvsimu', 'Hinvsimu');
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
solver = EITFEM(fm);
solver.sigmamin = 1e-9;
solver.zeta = 1e-9*ones(nel,1);
solver.mode = 'potential';
solver.vincl = vincl;
solver.Iel = Imeas_i;
solver.Uel = Umeas_i;
solver.SetInvGamma(1e-3, 3e-2);

solver_nld = EITFEM_NLD(fm);
solver_nld.sigmamin = 1e-9;
solver_nld.zeta = 1e-9*ones(nel,1);
solver_nld.mode = 'potential';
solver_nld.vincl = vincl;
solver_nld.Iel = Imeas;
solver_nld.Uel = Umeas;
solver_nld.SetInvGamma(1e-3, 3e-2);
disp('Forward problem solver set up')

%Set up the Total Variation prior:
%Start by calculating the value for alpha based on Gerardos method
ec = 10;%expected change
TVPrior = PriorTotalVariation(ginv, Hinv, [ec; ec]);
TVPrior.sigmaind = [1; 2];
disp('TV prior set up');

%Set up the positivity prior:
PosiPrior = PriorPositivityParabolic([1e-5; 1], [100; 0]);
disp('Positivity prior set up');

%Set up the plotter:
plotter = Plotter(ginv, Hinv);
plotter.plottypes = [1; 1];

%Finally, set up the inverse problem solver, in this case GN with
%linesearch:
resobj = cell(4, 1);
resobj{1} = solver;
resobj{2} = solver_nld;
resobj{3} = TVPrior;
resobj{4} = PosiPrior;
InvSolver = SolverLinesearch(resobj);
InvSolver.maxIter = 100;
InvSolver.Plotter = plotter;

%Make the initial guess and start solving!
sigmainitial = ones(size(ginv,1),1);
sigmainitial = Estimate_vec({sigmainitial; 0*sigmainitial});
disp('All set! Beginning to solve the inverse problem.')
reco = InvSolver.solve(sigmainitial);


