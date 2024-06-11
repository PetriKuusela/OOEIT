close all; clear all; clc
%This is the main-script for a simple traditional watertank EIT
%simulation WITH inverse crime demonstrating using the second order mesh.

InitializeOOEIT();%Make sure all the relevant paths are added

%% Create simulated data
disp('Creating simulated data');

simuname = '2ndsimu';
nel = 16;
vincl = true(nel,nel);%this variable tells which measurements are to be used; now its all of the measurements
vincl = vincl(:);

load_data = 0;%Do we run the simulation or just load the data.
if load_data
    load([simuname 'data.mat']);
    disp(['Loaded simulated data from file ' simuname 'data.mat']);
else
    load('meshfiles/watertank_mesh_2nd.mat');%Load mesh variables
    ginvsimu = ginv(:,1:2);%The node coordinates of the inverse mesh (a 1st order mesh where the conductivity is defined)
    Hinvsimu = Hinv;%The element connectivity of the inverse mesh
    gsimu = g;%The node coordinates of the forward mesh (a 2nd order mesh)
    Hsimu = H;%The element connectivity of the forward mesh
    elfacessimu = elfaces;%boundary element connectivity for each electrode
    fmsimu = ForwardMesh2nd(gsimu, Hsimu, elfacessimu);%create the mesh object
    imesh.g = ginvsimu;%also create an inverse mesh struct
    imesh.H = Hinvsimu;
    fmsimu.SetInverseMesh(imesh);%and pass the inverse mesh to the forward mesh object

    sigmasimu = GenerateEllipse(ginvsimu, 1, 10, 4e-2, 3e-2, 5e-2, 0e-2, 1e-2);%Generate a single blob as a target
    [Umeas, Imeas] = Simulation(fmsimu, [], sigmasimu, [], [], [], [0 0 0 0]);
    save([simuname 'data.mat'], 'Umeas', 'Imeas', 'sigmasimu', 'ginvsimu', 'Hinvsimu');
    disp(['simulation run succesfully, results saved in ' simuname 'data.mat']);
end

%% Load meshes for inverse problem
%Load meshes for solving the inverse problem, i.e. a 3D 2nd order forward mesh and a
%2D inverse mesh. (Here a same mesh is used as in simulations, i.e. inverse
%crime is committed)
load('meshfiles/watertank_mesh_2nd.mat');
ginv = ginv(:,1:2);
fm = ForwardMesh2nd(g, H, elfaces);
imesh.g = ginv;
imesh.H = Hinv;
fm.SetInverseMesh(imesh);
disp('meshfiles loaded');


%% Setting up the inverse solver

%Set up the forward problem solver:
solver = EITFEM(fm);
solver.mIncl = vincl;
solver.Iel = Imeas;
solver.Uel = Umeas;
solver.SetInvGamma(1e-3, 3e-2);
disp('Forward problem solver set up')

%Set up the Total Variation prior:
ec = 10;%expected change
TVPrior = PriorTotalVariation(ginv, Hinv, ec);
disp('TV prior set up');

%Set up the smoothness prior:%Uncomment and give this object to InvSolver if you want to use this
%SmoothPrior = PriorSmoothness(ginv, 0.05, 1, ones(size(ginv,1),1));
%disp('Smoothness prior set up');

%Set up the positivity prior:
PosiPrior = PriorPositivityParabolic(1e-5, 100);
disp('Positivity prior set up');

%Set up the plotter:
plotter = Plotter(ginv, Hinv);

%Finally, set up the inverse problem solver, in this case GN with
%linesearch:
resobj = cell(3, 1);
resobj{1} = solver;
resobj{2} = TVPrior;
resobj{3} = PosiPrior;
invSolver = SolverGN(resobj);
invSolver.maxIter = 100;
invSolver.plotter = plotter;

%Make the initial guess and start solving!
sigmainitial = ones(size(ginv,1),1);
disp('All set! Beginning to solve the inverse problem.')
reco = invSolver.Solve(sigmainitial);


