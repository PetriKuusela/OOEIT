classdef SolverLinesearch < handle
%SolverLinesearch class contains a Gauss-Newton optimization algorithm
%designed to be used in solving the inverse problem of EIT. It should
%function as well with any other optimization problem, with at most minor
%changes required. SolverLinesearch.solve(sigma) tries to find the minimum
%of the sum of the optimization functions set in the property ofuns. All
%the elements of ofuns have to have the methods OptimizationFunction(sigma)
%and GetHessAndGrad(sigma), the former giving the value that is to be
%minimized, and the latter giving the Hess-matrix and gradient of that
%value.
%
%Further, the first element of ofuns can have a plot(sigma) method to plot
%the fitting of the data. In order to plot the estimate, a plotter
%class can be added (in the property "Plotter"), which has to have a
%function plot(sigma), which then plots the current estimate.
%
%Author: Petri Kuusela, 5.4.2024
    
    properties
        erel            %The relative maximum difference in residuals, used to stop each linesearch
        estep           %The minimum steplength, used to stop the whole algorithm.
        estop           %The relative maximum difference in residuals, used to stop the whole algorithm
        nStepsBack      %Number of steps back we are comparing to stop the whole algorithm
        ofuns           %The optimization functions to be minimized, have to have functions OptimizationFunction(sigma), and CalculateHandGrad(sigma)
        plotLinesearch  %A flag: do we want to plot the optimization function values during linesearch
        hLinesearch     %The figure handle for plotting the above
        plotIterations  %A flag: do we want to plot the estimate each iteration
        plotData        %A flag: do we want to plot the data fitting each iteration
        hIterations     %The figure handle for plotting the above (sorry about the mixed naming)
        plotConvergence %A flag: do we want to plot convergence plots
        cws             %Convergence plots window size
        hConvergence    %The figure handle for plotting the above
        plotUQ          %A flag: do we want to plot credible interval widths
        nstd            %how many standard deviations should the credible interval be?
        maxIter         %Maximum number of GN iterations
        maxIterInLine   %Maximum number of iterations in one direction, i.e. inside one linesearch
        Plotter         %A plotter object, which can plot the estimate via function plot(sigma)
        tracker         %A cell array of vectors to track the evolution of the optimization function and other tracked values at GN iterations
        showSplitVals   %A flag: plot each optimization function value separately during linesearch
        outputFunVals   %A flag: output the optimization function values on the console each iteration
        alwaysMove      %A flag: Do we move a little bit even if ofuns value does not decrease
        steepestDesc    %A flag: Use steepest descend instead of GN
        trueEst         %the true value of estimates, that can be known e.g. in simulations
    end
    
    methods
        function obj = SolverLinesearch(ofuns)
            %Class contstructor. Input is the cell array containing
            %optimization function -objects.
            obj.ofuns = ofuns;
            obj.erel = 1e-4;%The relative maximum difference in residuals, used to stop the linesearch
            obj.estep = 1e-7;%The maximum steplength, used to stop the whole algorithm.
            obj.estop = 1e-2;
            obj.nStepsBack = 20;
            obj.plotLinesearch = 1;
            obj.plotIterations = 1;
            obj.plotData = 1;
            obj.plotConvergence = 1;
            obj.maxIter = 10;
            obj.maxIterInLine = 30;
            obj.showSplitVals = 0;
            obj.outputFunVals = 1;
            obj.alwaysMove = 1;
            obj.cws = -2;
            obj.steepestDesc = 0;
            obj.nstd = 2;
        end
        
        function reco = solve(self, sigest)
            %The main function of the object, this is where the
            %Gauss-Newton is implemented. The argument "sigest" is the
            %initial guess for the estimate.
            
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
                if isempty(self.hIterations)
                    self.hIterations = figure();
                end
                self.ofuns{1}.Plot(sigest)
            end
            if self.plotConvergence && isempty(self.hConvergence)
                self.hConvergence = figure();
            end
            self.tracker = cell(length(self.ofuns)+1,1);
                
            
            iter = 1; %The number of current iteration
            self.UpdateTracker(sigest, iter, 0);%Track the evolution of the optimization functions
            laststep = 1;%adjust the initial guess for step length based on where the last minimum was found
            
            while (iter <= self.maxIter)
                %This is the main loop for GN
                
                %Calculate the gradients and Hessians
                nobj = length(self.ofuns);
                Hess = 0;
                grad = 0;
                for io = 1:nobj
                    [tempHess, tempgrad] = self.ofuns{io}.GetHessAndGrad(sigest);
                    Hess = Hess + tempHess;
                    grad = grad + tempgrad;
                end
                if self.plotUQ
                    if isa(Hess, 'Estimate_Hess')
                        invHess = inv(Hess.Matrix());
                    else
                        invHess = inv(Hess);
                    end
                    self.Plotter.plotUQ(self.nstd*sqrt(diag(invHess)));
                    self.Plotter.lineplot(sigest, self.nstd*sqrt(diag(invHess)), self.trueEst);
                end

                %Line search (assume there is only one local minimum on the region
                %accessible to this algorithm)

                fval = self.CalculateFval(sigest); %the optimization function value of current point
                if self.steepestDesc
                    deltasigma = -grad; %Search direction of steepest descend, i.e. negative gradient.
                else
                    deltasigma = -Hess\grad; %The direction determined by the Gauss Newton algorithm
                end
                if self.showSplitVals
                    fvalsplit = self.SplitFvals(sigest)';
                end

                %In this linesearch we compare four points at a time in order to determine where
                %to try the next point; the leftmost, rightmost and the ceter points,
                %and in addition the new lastly calculated point which is between the
                %center and either the left or right. From these points the left or
                %right point will be left out and the remaining points will be renamed
                %according to their position, and then a new point will be tested.

                cont = 1;%continue flag for the linesearch
                d = [0, laststep]';%vector containing different relative step lengths tried
                ii = 2;%the line search iteration counter
                r = 0; %the index of the right point (0 = we're still checking if we want to go even further right)
                l = 1; %the index of the left point
                c = 1; %the index of the center point
                while cont == 1%This is the loop for the linesearch
                    sigest_new = sigest + d(ii)*deltasigma;
                    fval(ii) = self.CalculateFval(sigest_new);%in each point we calculate the optimization function value
                    if self.showSplitVals
                        fvalsplit(ii,:) = self.SplitFvals(sigest_new)';
                    end
                    if r == 0 %we're still going right
                        if fval(ii) < fval(ii-1) %continue right
                            l = ii-1;
                            d(ii+1) = 2*d(ii);%The new step is double the last one
                        else%stop going right
                            r = ii;
                            c = ii-1;
                            d(ii+1) = 0.5*(d(r)+d(c));%the new step is the avg of the last two tried
                        end
                    else%the minimum is between our left and right points
                        if (d(ii)<d(c) && fval(ii) < fval(c))%The minimum is left from the point c
                            r = c;
                            c = ii;
                            d(ii+1) = 0.5*(d(c)+d(r));
                        elseif (d(ii)<d(c) && fval(ii) > fval(c))%The minimum might be right from the point c
                            l = ii;
                            d(ii+1) = 0.5*(d(c)+d(r));
                        elseif (d(ii)>d(c) && fval(ii) < fval(c))%The minimum is right from the point c
                            l = c;
                            c = ii;
                            d(ii+1) = 0.5*(d(c)+d(r));
                        elseif (d(ii)>d(c) && fval(ii) > fval(c))%The minimum might be left from the point c
                            r = ii;
                            if (c~=l)
                                d(ii+1) = 0.5*(d(l)+d(c));
                            else%c==1 means we have not found a point with lower fval than the start point d(1) = 0.
                                d(ii+1) = 0.5*(d(l)+d(ii));
                            end
                        else%we should not end up here
                            error('Something wrong with line search');
                        end
                    end
                    
                    %Plot the fvals of the linesearch
                    if self.plotLinesearch
                        if self.showSplitVals
                            self.LinesearchPlot(d, fvalsplit, ii);
                        else
                            self.LinesearchPlot(d, fval, ii);
                        end
                    end
                    if (r > 0)%see if we're happy with our results and can jump out of the linesearch
                        if (abs(fval(l) - min(fval)) < self.erel*fval(ii) && abs(fval(r) - min(fval)) < self.erel*fval(ii))
                            cont = 0;
                        end
                    end
                    ii = ii +  1;
                    if (ii > self.maxIterInLine)
                        %Reached maximum number of iterations in
                        %linesearch, so moving on.
                        warning('Line search did not converge, continue anyways.');%The values did not converge in the maximum number of iterations
                        cont = 0;
                    end
                    
                end%end linesearch, move on to the next GN iteration
                
                if self.alwaysMove
                    if (min(fval(2:end)) > fval(1))%we want to move always, so not choosing the initial point
                        warning(['The residual increased on iteration ', num2str(iter)]);%choose the best point of the ones tried
                    end

                    [minfval, minii] = min(fval(2:end));%find the index of the minimum. (we use minii+1 as the index below because minii is the index in fval(2:end))
                    sigest = sigest + d(minii+1)*deltasigma;%update the conductivity distribution
                    laststep = 1;% d(minii+1);%Save the last steplength to start with that as a guess in the next iteration
                else%We have the option to not move at all (this means stopping the GN)
                    [minfval, minii] = min(fval);%find the index of the minimum.
                    sigest = sigest + d(minii)*deltasigma;%update the conductivity distribution
                    laststep = 1;% d(minii+1);%Save the last steplength to start with that as a guess in the next iteration
                end

                if self.plotIterations%Plot the estimate
                    self.Plotter.plot(sigest)
                end
                if self.plotData%Plot how the data fits
                    set(0, 'CurrentFigure', self.hIterations);self.ofuns{1}.Plot(sigest)%plot how the measurements align with the simulated results for current estimation
                    self.hIterations.Name = ['Iteration ', num2str(iter)];%mark down the iteration too on this window
                end
                if self.outputFunVals
                    disp(['Finished iteration ' num2str(iter) ' with optimization function value ' num2str(minfval, '%e')]);
                end

                iter = iter + 1;
                self.UpdateTracker(sigest, iter, norm(d(minii)*deltasigma));

                if self.plotConvergence
                    self.ConvergencePlots(self.cws);
                end
                drawnow;
                
                if (norm(d(minii+1)*deltasigma) < self.estep)
                    disp('moving stopped with the step length of:');
                    disp(norm(d(minii+1)*deltasigma));
                    break;
                end

                if iter > self.nStepsBack && min(fval(2:end)) > (1-self.estop)*self.tracker{length(self.ofuns)+1}(iter-self.nStepsBack)
                   disp('End Gauss-Newton, residual converged');
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
            res = 0;
            for io = 1:length(self.ofuns)
                tempres = self.ofuns{io}.OptimizationFunction(sigest);
                self.tracker{io}(iter) = tempres;
                res = res + tempres;
            end
            self.tracker{length(self.ofuns)+1}(iter) = res;
            self.tracker{length(self.ofuns)+2}(iter) = dsnorm;
        end

        function LinesearchPlot(self, d, fval, ii)
            set(0, 'CurrentFigure', self.hLinesearch);
            if self.showSplitVals
                legends = cell(length(self.ofuns)*2,1);
                for io = 1:length(self.ofuns)
                    legends{2*io-1} = num2str(io);
                    legends{2*io} = '';
                end
                [dsorted, sorti] = sort(d(1:end-1));
                fvalsorted = fval(sorti,:);
                clf;
                hold on;
                for io = 1:length(self.ofuns)
                    plot(dsorted, fvalsorted(:,io), 'o-');
                    plot(d(ii), fval(ii,io), 'r+');
                end
                hold off;
                legend(legends);
            else
                plot(d(1:end-1), fval, 'bo', d(ii), fval(ii), 'r+'); %plot the progress of the linesearch
            end
            self.hLinesearch.Name = ['Linesearch iteration ', num2str(ii)];%mark down the iteration too on this window
            drawnow;
        end

        function ConvergencePlots(self, windowsize)
            set(0, 'CurrentFigure', self.hConvergence);
            clf;
            if length(self.ofuns) < 3
                tiledlayout(2,2);
            else
                tiledlayout(ceil((length(self.ofuns)+2)/3), 3);
            end

            for io = 1:length(self.tracker)
                nexttile;
                if windowsize == 0
                    plot(0:length(self.tracker{io})-1, self.tracker{io});
                elseif windowsize < 0
                    plot(min(length(self.tracker{io})-1,-windowsize):(length(self.tracker{io})-1), self.tracker{io}(min(-windowsize+1,length(self.tracker{io})):end));
                else
                    plot(max(0,length(self.tracker{io})-windowsize):(length(self.tracker{io})-1), self.tracker{io}(max(1,length(self.tracker{io})-windowsize+1):end));
                end
                if io <= length(self.ofuns)
                    title(class(self.ofuns{io}));
                elseif io == length(self.ofuns)+1
                    title('Optimization function');
                elseif io == length(self.ofuns)+2
                    title('Reconstruction change');
                end
            end

        end
        
        function [res, grads, Hs] = SplitFvals(self, sigest)
            %The first output is the same as CalculateFval above, but
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
        
        
    end    
    
    
end