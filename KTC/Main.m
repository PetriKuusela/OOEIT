close all; clear all; clc;

addpath('..');
initialize_OOEIT();

%% Load meshes for inverse problem
%Load meshes for solving the inverse problem, i.e. a 3D forward mesh and a
%2D inverse mesh.
load('meshfiles/Mesh_sparse.mat');
%ginv = Mesh.g;
g = Mesh.g;
H = Mesh.H;
elfaces = Mesh.elfaces;
nel = length(elfaces);
fm = ForwardMesh1st(g, H, elfaces);
%imesh.g = ginv;
%imesh.H = Hinv;
%fm.SetInverseMesh(imesh);
disp('meshfiles loaded');


%% Setting up the inverse solver

%Set up the forward problem solver:
solver = EITFEM(fm);
solver.sigmamin = 1e-9;
solver.zeta = 1e-9*ones(nel,1);
solver.mode = 'current';
disp('Forward problem solver set up')

%Set up the Total Variation prior:
%Start by calculating the value for alpha based on Gerardos method
ec = 1;%expected change
TVPrior = PriorTotalVariation(g, H, ec);
disp('TV prior set up');

%Set up the smoothness prior:
SmoothPrior = PriorSmoothness(g, 0.2, 0.1, ones(size(g,1),1));
disp('Smoothness prior set up');


%Set up the white noise prior:
%WNPrior = PriorWhiteNoise(0.7, 0.5);

%Set up the positivity prior:
%PosiPrior = PriorPositivityParabolic([1e-5; 1e-8], [100, 1]);
PosiPrior = PriorPositivityParabolic(1e-5, 100);
disp('Positivity prior set up');

%Set up the plotter:
plotter = Plotter(g, H);

%Finally, set up the inverse problem solver, in this case GN with
%linesearch:
resobj = cell(3, 1);
resobj{1} = solver;
resobj{2} = TVPrior;
resobj{3} = PosiPrior;
InvSolver = SolverGN(resobj);
InvSolver.maxIter = 150;
InvSolver.Plotter = plotter;

sigmainit = ones(size(g,1),1);
hesti = Estimate_vec({ones(size(g,1),1), 1e-5*ones(32,1)});
InvMpat = tril(-ones(nel),-1);
InvMpat(:,end) = [];
score = zeros(7, 3);
for ilevel = 1:7
    refname = ['EvaluationData/evaluation_datasets/level' num2str(ilevel) '/ref.mat'];
    load(refname);

    solver.eps = 0;
    solver.Iel = Injref(:);
    usempat = 1;
    if usempat == 0
        Uelrefmat = reshape(Uelref, 31, []);
        Uelrefvec = InvMpat*Uelrefmat;
        Uelrefvec = Uelrefvec - mean(Uelrefvec, 'omitnan');
        Uelrefvec = Uelrefvec(:);
        solver.vincl = ~isnan(Uelrefvec);
        Uelrefvec(isnan(Uelrefvec)) = [];
        solver.Uel = Uelrefvec;
    else
        solver.vincl = ~isnan(Uelref);
        solver.Uel = Uelref(~isnan(Uelref));
        solver.Mpat = Mpat;
    end
    solver.SetInvGamma(1e-3, 3e-2);
    solver.zind = 2;
    [eps, hest] = SolveEpsilonCorrection_withz(solver, hesti);
    solver.zeta = hest.estimates{2};

    for idata = 1:3

        dataname = ['EvaluationData/evaluation_datasets/level' num2str(ilevel) '/data' num2str(idata) '.mat'];
        load(dataname);
        
        solver.Iel = Inj(:);
        Uelmat = reshape(Uel, 31, []);
        if usempat == 0
            Uelvec = InvMpat*Uelmat;
            Uelvec = Uelvec - mean(Uelvec, 'omitnan');
            Uelvec = Uelvec(:);
            solver.vincl = ~isnan(Uelvec);
            eps2 = eps(solver.vincl);
            Uelvec(isnan(Uelvec)) = [];
            solver.Uel = Uelvec(:);
        else
            solver.vincl = ~isnan(Uel);
            eps2 = eps(solver.vincl);
            solver.Uel = Uel(~isnan(Uel));
            solver.eps = eps2;
        end
        solver.SetInvGamma(1e-3, 3e-2);

        disp('All set! Beginning to solve the inverse problem.')
        reco = InvSolver.solve(hest.estimates{1}*sigmainit);

        save(['Figures/data' num2str(ilevel) '_' num2str(idata) '.mat'], 'reco', 'InvSolver');
        %saveas(InvSolver.Plotter.fig, ['Figures/figure' num2str(ilevel) '_' num2str(idata) '.png']);

        pixreco = interpolateRecoToPixGrid(reco, fm);
        tf = figure(); imagesc(pixreco);
        saveas(tf, ['Figures/figure' num2str(ilevel) '_' num2str(idata) '.png']);

        [level, x] = Otsu2(pixreco, 256);

        pixreco(pixreco<x(level(1))) = 0;
        pixreco(pixreco<x(level(2)) & pixreco > 0) = 0.5e-5;
        pixreco(pixreco > 1e-5) = 1e-5;
        pixreco = pixreco*1e5;

        numpix0 = sum(pixreco <0.25, 'all');
        numpix05 = sum(pixreco>0.25 & pixreco<0.75, 'all');
        numpix1 = sum(pixreco>0.75, 'all');

        if numpix0 > numpix05 + numpix1
            pixreco(pixreco > 0.25) = 1;
            pixreco(pixreco < 0.25) = 0.5;
        elseif numpix1 > numpix05 + numpix0
            pixreco(pixreco < 0.75) = 0;
            pixreco(pixreco > 0.75) = 0.5;
        end

        pixreco(pixreco < 0.25) = 10;
        pixreco(pixreco < 0.75) = 0;
        pixreco(pixreco > 0.75 & pixreco < 1.25) = 2;
        pixreco(pixreco > 5) = 1;

        imagesc(pixreco);
        saveas(tf, ['Figures/figure' num2str(ilevel) '_' num2str(idata) '_seg.png']);
        close(tf);

        load(['EvaluationData/GroundTruths/level_' num2str(ilevel) '/' num2str(idata) '_true.mat']);
        score(ilevel, idata) = scoringFunction(truth, pixreco);

    end
end

save('Figures/scores.mat', 'score');
