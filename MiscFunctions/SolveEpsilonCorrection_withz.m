function [eps, hest] = SolveEpsilonCorrection_withz(solver, sigmainit)

    if nargin < 2
        sigmainit = Estimate_vec({1; 1e-5*ones(solver.fmesh.nel,1)});
    end
    
    hw = HomogeneousWrapper(solver);
    PosiPrior = PriorPositivityParabolic([1e-5; 1e-8], [1e2; 1]);
    resobj = cell(2,1);
    resobj{1} = hw;
    resobj{2} = PosiPrior;
    InvSolver = SolverGN(resobj);
    InvSolver.maxIter = 20;
    InvSolver.plotIterations = 1;
    InvSolver.showSplitVals = 1;
    
    hinit = sigmainit;
    hinit.estimates{1} = 1;
    hest = InvSolver.solve(hinit);
    
    if strcmp(solver.mode, 'potential')
        Iel = hw.SolveForwardVec(hest);
        eps = solver.Iel-Iel;
    elseif strcmp(solver.mode, 'current')
        Uel = hw.SolveForwardVec(hest);
        eps = solver.Uel-Uel;
    else
        error('Unrecognized solver mode');
    end

end