classdef SolverGN < handle
%SolverGN class contains a Gauss-Newton optimization algorithm
%designed to be used in solving the inverse problem of EIT. It should
%function as well with any other optimization problem, with at most minor
%changes required. SolverGN.solve(sigma) tries to find the minimum
%of the sum of the optimization functions set in the property ofuns. All
%the elements of ofuns have to have the methods OptimizationFunction(sigma)
%and GetHessAndGrad(sigma), the former giving the value that is to be
%minimized, and the latter giving the Hess-matrix and gradient of that
%value.
%
%Further, the ofuns objects can have plot(sigma) methods to plot
%the fitting of the data. In order to plot the estimate, a plotter
%class can be added (in the property "Plotter"), which has to have a
%function plot(sigma), which then plots the current estimate.
%
%Author: Petri Kuusela, 31.5.2024
    
    properties
        erel            %The relative maximum difference in residuals, used to stop each linesearch
        estep           %The minimum steplength advance in the GN iteration, used to stop the whole algorithm if it is not progressing anywhere
        estop           %The relative maximum difference in residuals compared to nStepsBack steps back, used to stop the whole algorithm
        nStepsBack      %Number of steps back we are comparing to stop the whole algorithm
        ofuns           %The optimization functions to be minimized, have to have functions OptimizationFunction(sigma), and GetHessAndGrad(sigma)
        plotLinesearch  %A flag: do we want to plot the optimization function values during linesearch
        hLinesearch     %The figure handle for plotting the above
        plotIterations  %A flag: do we want to plot the estimate each iteration
        plotData        %A flag: do we want to plot the data fitting each iteration
        plotDataInd     %The indices of ofuns that we want to call the plot function
        hData           %The figure handle for plotting the above
        plotConvergence %A flag: do we want to plot convergence plots
        cws             %Convergence plots iteration-window size (negative value means show all values except -cws from beginning. positive values mean show cws last values. 0 means show all)
        hConvergence    %The figure handle for plotting the above
        plotUQ          %A flag: do we want to plot credible interval widths
        nstd            %how many standard deviations should the credible interval be?
        maxIter         %Maximum number of GN iterations
        maxIterInLine   %Maximum number of iterations inside each linesearch
        Plotter         %A plotter object, which can plot the estimate via function Plotter.plot(sigma)
        tracker         %A cell array of vectors to track the evolution of the optimization function and other tracked values at GN iterations
        showSplitVals   %A flag: plot each optimization function value separately during linesearch
        outputFunVals   %A flag: output the optimization function values on the console each iteration
        alwaysMove      %A flag: Do we move a little bit even if ofuns value does not decrease
        steepestDesc    %A flag: Use steepest descend for search direction instead of GN. May be used in extremely non-linear locations, where GN fails.
        parabolic       %A flag: Try to find the minimum in linesearch using parabolic fits. This usually speeds up the computations, but may cause problems in some cases
        trueEst         %the true value of estimates, that can be known e.g. in simulations. Can be used by the plotter for lineplots.
        uselastd        %A flag, start the linesearch at relative distance of last minimum. Sometimes speeds up computations, sometimes causes issues.
    end
    
    methods
        function obj = SolverGN(ofuns)
            %Class contstructor. Input is a cell array containing
            %optimization function -objects.
            obj.ofuns = ofuns;

            %For explanations of the properties, please check the
            %properties-section. Populate with default values:
            obj.erel = 1e-5;
            obj.estep = 1e-8;
            obj.estop = 1e-2;
            obj.nStepsBack = 10;
            obj.plotLinesearch = 1;
            obj.plotIterations = 1;
            obj.plotData = 1;
            obj.plotDataInd = 1;
            obj.plotConvergence = 1;
            obj.maxIter = 100;
            obj.maxIterInLine = 30;
            obj.showSplitVals = 0;
            obj.outputFunVals = 1;
            obj.alwaysMove = 1;
            obj.cws = -2;
            obj.steepestDesc = 0;
            obj.nstd = 2;
            obj.parabolic = 1;
            obj.uselastd = 1;
        end
        
        function reco = solve(self, sigest)
            %The main function of the object, this is where the
            %Gauss-Newton is implemented. The argument "sigest" is the
            %initial guess for the estimate, and output "reco" is the final
            %reconstruction obtained.
            
            %First check that all the relevant figures exists, and if not,
            %create them:
            if self.plotLinesearch && isempty(self.hLinesearch)
                self.hLinesearch = figure();%create a new figure if one does not exist
            end
            if self.plotIterations || self.plotUQ
                if isempty(self.Plotter)
                    warning('No Plotter set, so could not plot the iterations nor the UQ');
                    self.plotIterations = 0;
                    self.plotUQ = 0;
                end
            end
            if self.plotData
                if isempty(self.hData)
                    self.hData = figure();
                end
                self.ofuns{1}.Plot(sigest)%Plot the initial data fit
            end
            if self.plotConvergence && isempty(self.hConvergence)
                self.hConvergence = figure();
            end
            self.tracker = cell(length(self.ofuns)+1,1);%Store the initial values for the tracker
                
            
            iter = 1; %The number of current iteration
            lastd = 1; %lastd contains the relative distance of minimum in last linesearch, if self.uselastd == 1. Start with lastd = 1
            self.UpdateTracker(sigest, iter, 0);%Track the evolution of the optimization functions
            nobj = length(self.ofuns);
            
            while (iter <= self.maxIter)
                %This is the main loop for GN
                
                %Calculate the gradients and Hessians
                Hess = 0;
                grad = 0;
                for io = 1:nobj %The final gradient and Hessian are just sums of the gradient and Hessian of each optimization function
                    [tempHess, tempgrad] = self.ofuns{io}.GetHessAndGrad(sigest);
                    Hess = Hess + tempHess;
                    grad = grad + tempgrad;
                end

                
                if self.plotUQ %Plot the credible interval widths and given lineplots
                    %Using Gaussian approximation, the credible interval
                    %widths depend on the inverse of the Hessian of the
                    %posterior distribution around the MAP.
                    if isa(Hess, 'Estimate_Hess')
                        invHess = inv(Hess.Matrix());
                    else
                        invHess = inv(Hess);
                    end
                    self.Plotter.plotUQ(self.nstd*sqrt(diag(invHess)));
                    %The following lineplot plots the estimate and credible
                    %interval on a line through the imaging domain. Check
                    %the function comments for required changes to
                    %configure it matching your needs:
                    self.Plotter.lineplot(sigest, self.nstd*sqrt(diag(invHess)), self.trueEst);
                end

                %Compute the search direction:
                if self.steepestDesc
                    deltasigma = -grad; %Search direction of steepest descend, i.e. negative gradient.
                else
                    deltasigma = -Hess\grad; %The direction determined by the Gauss Newton algorithm
                end

                %If the optimizaiton function values are plotted
                %separately, collect the necessary values at start point
                if self.showSplitVals
                    fvalsplit = self.SplitFvals(sigest)';
                end

                %Line search begins (assume there is only one local minimum on the region
                %accessible to this algorithm)

                %This class contains two possible options for the
                %linesearch. The default choice is parabolic fitting, where
                %a parabola is fitted to the 5 lowest ofun value points,
                %and the next point to test is at the minimum of this
                %parabola. If this yields bad results, the algorithm
                %automatically reverts to non-parabolic search, where a
                %location between the current minimum point and the higher
                %one of its neighbours is tested next.

                %if self.parabolic == 1, then each linesearch starts with
                %parabolic search, even if it was abandoned on last search.
                
                parabolicnow = self.parabolic;
                notgood = 0;%This flag determines if we need to abandon parabolic search
                cont = 1;%continue flag for the linesearch
                d = [-lastd 0 lastd]';%vector containing different relative step lengths tried
                ii = 3;%the line search iteration counter (the first three points are fixed)
                fval = zeros(3,1);
                fval(2) = self.CalculateFval(sigest); %the optimization function value of the fixed points
                fval(1) = self.CalculateFval(sigest + d(1)*deltasigma);
                fval(3) = self.CalculateFval(sigest + d(3)*deltasigma);
                while cont == 1%This is the loop for the linesearch
                    
                    %find the minimum point
                    [~,minind] = min(fval);
                    if minind == 1%Check if our minimum might be on the left of any of the current points
                        direction = -1;
                    elseif minind == ii%Check if our minimum might be on the right of any of the current points
                        direction = 1;
                    else%The minimum is somewhere between the leftmost and rightmost point
                        direction = 0;
                    end
                    dists = abs(d-d(minind));
                    [~,inds] = sort(dists);
                    inds = inds(1:min(length(inds),5));%These are the indices of 5 lowest ofun values
                    maxd = dists(inds(end));
                    inds2 = max(minind-2,1):min(minind+2,ii);%These are the indices of the 5 point around the current minimum index
                    if parabolicnow%Do the parabolic fit to find next point to try
                        [d(ii+1), notgood] = self.FindParabolicMinimum(d(inds), fval(inds));
                    end
                    if ~parabolicnow || notgood %No parabolic fit
                        if direction == -1%Our minimum is somewhere to the left
                            d(ii+1) = 2*d(1);
                        elseif direction == 1%Our minimum is somewhere to the right
                            d(ii+1) = 2*d(ii);
                        else%Our minimum is between our left- and rightmost points
                            %We test next the point halfway between current
                            %minimum and the one of its neighboring points,
                            %whose ofun value is higher
                            if fval(minind-1) > fval(minind+1)
                                d(ii+1) = 0.5*(d(minind-1)+d(minind));
                            else
                                d(ii+1) = 0.5*(d(minind+1)+d(minind));
                            end
                        end
                    end

                    ii = ii + 1;%iteration counter

                    %The new point to try
                    sigest_new = sigest + d(ii)*deltasigma;
                    fval(ii) = self.CalculateFval(sigest_new);%in each point we calculate the optimization function value
                    if self.showSplitVals
                        fvalsplit(ii,:) = self.SplitFvals(sigest_new)';%Get also the separate optimization function values of each function
                    end
                    if parabolicnow && ((~direction && fval(ii) > max(fval(inds2))) || notgood || abs(d(ii)-d(minind))>maxd && fval(ii)>fval(minind))
                        %We have a reason to abandon the parabolic-fitting
                        %linesearch (will be reset at next linesearch again)
                        parabolicnow = 0;
                    end

                    %Plot the fvals of the linesearch
                    if self.plotLinesearch
                        if self.showSplitVals
                            self.LinesearchPlot(d, fvalsplit, ii);
                        else
                            self.LinesearchPlot(d, fval, ii);
                        end
                    end
                    
                    %see if we're happy with our results and can jump out of the linesearch
                    if (mean(fval(inds) - min(fval)) < self.erel*min(fval))
                        cont = 0;
                    end
                    if (ii > self.maxIterInLine)
                        %Reached maximum number of iterations in
                        %linesearch, so moving on.
                        %If you want to monitor how well the linesearch
                        %converges, you can comment the next line to get a
                        %warning each time it did not converge:
                        %warning('Line search did not converge, continue anyways.');%The values did not converge in the maximum number of iterations
                        cont = 0;
                    end

                    %sort our distances (only the last value should be out
                    %of its place
                    [d, dinds] = sort(d);
                    %the fvals have to match the distance indices as well
                    fval = fval(dinds);
                    
                end%end linesearch, move on to the next GN iteration
                

                if self.alwaysMove
                    preval = fval(d==0);%the option d = 0 is not permitted
                    fval(d==0) = [];
                    d(d==0) = [];
                    if (min(fval) > preval)%we want to move always, so not choosing the initial point
                        warning(['The residual increased on iteration ', num2str(iter)]);%choose the best point of the ones tried
                    end

                    [minfval, minii] = min(fval);%find the index of the minimum.
                    sigest = sigest + d(minii)*deltasigma;%update the conductivity distribution
                else%We have the option to not move at all (this means stopping the GN)
                    [minfval, minii] = min(fval);%find the index of the minimum.
                    sigest = sigest + d(minii)*deltasigma;%update the conductivity distribution
                end
                if self.uselastd%we want to start the next linesearch at the relative distance of the minimum of this linesearch
                    lastd = d(minii);
                else
                    lastd = 1;
                end

                if self.plotIterations%Plot the estimate
                    self.Plotter.plot(sigest)
                end
                if self.plotData%Plot how the data fits
                    for iplot = 1:length(self.plotDataInd)
                        if length(self.hData) < iplot || isempty(self.hData(iplot))
                            self.hData(iplot) = figure();
                        end
                        set(0, 'CurrentFigure', self.hData);self.ofuns{self.plotDataInd(iplot)}.Plot(sigest)%plot how the measurements align with the simulated results for current estimation
                        self.hData.Name = ['Iteration ', num2str(iter)];%mark down the iteration too on this window
                        title(['Data fit, iteration: ' num2str(iter)]);
                    end
                end
                if self.outputFunVals
                    disp(['Finished iteration ' num2str(iter) ' with optimization function value ' num2str(minfval, '%e')]);
                end

                %GN iteration counter
                iter = iter + 1;
                %Track ofun values for convergence plots
                self.UpdateTracker(sigest, iter, norm(d(minii)*deltasigma));

                if self.plotConvergence%Plot the convergence plots
                    self.ConvergencePlots(self.cws);
                end
                drawnow;
                
                if (norm(d(minii)*deltasigma) < self.estep*norm(sigest))%The GN is practically not moving
                    disp('moving stopped with the step length of:');
                    disp(norm(d(minii)*deltasigma));
                    break;
                end

                if iter > self.nStepsBack && min(fval) > (1-self.estop)*self.tracker{length(self.ofuns)+1}(iter-self.nStepsBack)
                   disp('End Gauss-Newton, residual converged');%The GN has converged properly
                   break;
                end

            end %next iteration of GN

            %output the reason for halting the algorithm if maximum number
            %of iterations was reached
            if iter > self.maxIter-1
                disp ('Reached the limiting number of iterations in Gauss-Newton solver');
            end
            
            reco = sigest;%This is the output of the function
            
        end
        
        function res = CalculateFval(self, sigest)
            %Calculate the sum of optimization function values at sigest
            nobj = length(self.ofuns);
            res = 0;
            for io = 1:nobj
                res = res + self.ofuns{io}.OptimizationFunction(sigest);
            end
        end

        function UpdateTracker(self, sigest, iter, dsnorm)
            %Calculate all the optimization function values to the tracker
            %cell array
            res = 0;
            for io = 1:length(self.ofuns)
                tempres = self.ofuns{io}.OptimizationFunction(sigest);
                self.tracker{io}(iter) = tempres;
                res = res + tempres;
            end
            self.tracker{length(self.ofuns)+1}(iter) = res;%Store also the sum of the optimization functions
            self.tracker{length(self.ofuns)+2}(iter) = dsnorm;%and the norm of the change compared to last iteration
        end

        function LinesearchPlot(self, d, fval, ii)
            %This function plots the optimization function values as a
            %function of the relative distance d at the
            %points tried by the linesearch
            set(0, 'CurrentFigure', self.hLinesearch);
            if self.showSplitVals%Do we want to plot the different functionals separately
                %First make a cell array for the legend of the figure
                legends = cell(length(self.ofuns)*2,1);
                for io = 1:length(self.ofuns)
                    legends{2*io-1} = num2str(io);
                    legends{2*io} = '';%Every other plotted element will be just the mark on the current step of linesearch, thus no legend for that
                end

                %Sort the points by relative distance
                [dsorted, sorti] = sort(d(1:end-1));
                fvalsorted = fval(sorti,:);

                clf;
                hold on;
                for io = 1:length(self.ofuns)
                    plot(dsorted, fvalsorted(:,io), 'o-');%plot the values of a single optimization function as a line
                    plot(d(ii), fval(ii,io), 'r+');%and mark a red cross to the last tried point
                end
                hold off;
                legend(legends);
            else
                plot(d, fval, 'bo', d(ii), fval(ii), 'r+'); %plot the progress of the linesearch
            end
            self.hLinesearch.Name = ['Linesearch iteration ', num2str(ii)];%mark down the iteration too on this window
            title(['Linesearch, iteration:' num2str(ii)]);
            drawnow;
        end

        function ConvergencePlots(self, windowsize)
            %This function plots the values of different optimization
            %functions as a function of the GN iteration number. In
            %addition, the sum of the optimization functions and the change
            %compared to last iteration are plotted. These plots may be
            %useful in determining has the GN algorithm converged properly,
            %or to give hints on how should the parameter values be tuned.

            set(0, 'CurrentFigure', self.hConvergence);
            clf;

            %Try to automatically find a somewhat suitable tiling
            if length(self.ofuns) < 3
                t = tiledlayout(2,2);
            else
                t = tiledlayout(ceil((length(self.ofuns)+2)/3), 3);
            end

            %then start plotting
            for io = 1:length(self.tracker)
                nexttile;
                if windowsize == 0%In this case we want to plot all the iterations
                    plot(0:length(self.tracker{io})-1, self.tracker{io});
                elseif windowsize < 0%In this case we leave -windowsize of iterations from the beginning
                    plot(min(length(self.tracker{io})-1,-windowsize):(length(self.tracker{io})-1), self.tracker{io}(min(-windowsize+1,length(self.tracker{io})):end));
                else%In this case we plot the lase windowsize iterations
                    plot(max(0,length(self.tracker{io})-windowsize):(length(self.tracker{io})-1), self.tracker{io}(max(1,length(self.tracker{io})-windowsize+1):end));
                end

                %Make a title for the subfigure
                if io <= length(self.ofuns)%Take the title straight from the prior class name
                    title(class(self.ofuns{io}));
                elseif io == length(self.ofuns)+1%Plotting the sum of all functionals
                    title('Optimization function');
                elseif io == length(self.ofuns)+2%Plotting the changes between iterations
                    title('Reconstruction change');
                end
            end

            %Finally, make a title on the top
            title(t, 'Convergence plots');

        end
        
        function [res, grads, Hs] = SplitFvals(self, sigest)
            %The first output is the same as CalculateFval above, but this
            %returns a vector of same length as self.ofuns, with the fval
            %of each ofun separate as it's own element. The two further
            %outputs will not be calculated unless requested, and they give
            %the gradients and Hess-matrices of each ofun in cell arrays.
            nobj = length(self.ofuns);
            res = zeros(nobj,1);
            for io = 1:nobj
                res(io) = self.ofuns{io}.OptimizationFunction(sigest);
            end
            if nargout > 1%Gradients have been requested
                grads = cell(nobj,1);
                Hs = cell(nobj,1);
                for io = 1:nobj
                    [grads{io}, Hs{io}] = self.ofuns{io}.GetHessAndGrad(sigest);
                end
            end
        end
        
        function [minpos, notgood] = FindParabolicMinimum(self, x, y)
            %This function finds the minimum point of a parabola fit to data
            %[x, y] by LS. Further, the flag notgood is set to 1 if the
            %parabola opens downwards or the LS matrix is badly conditioned,
            %because then we know the parabolic
            %fit is useless for finding the minimum in linesearch.
            
            A = [x.^2 x ones(length(x),1)];
            ATA = A'*A;
            if rcond(ATA) < 1e-10%The parabolic fit will not be good
                notgood = 1;
                minpos = 0;%this will not be used in the linesearch
            else
                coefs = (ATA)\A'*y;%This is just a basic LS-fit
                if coefs(1) < 0
                    notgood = 1;%The parabola opens downwards, so it is not good for finding the minimum
                    minpos = 0;%this will not be used in the linesearch
                else
                    minpos = -0.5*coefs(2)/coefs(1);%Then find the analytic minimum of the fitted parabola
                    notgood = 0;
                end
            end
        end
        
    end    
    
    
end