classdef Plotter < handle
%A Plotter class that can be assigned to inverse problem solver (e.g.
%SolverGN.m) to plot the estimate in between the iterations.
    
    properties
        g       %The array of node co-ordinates
        H       %The elements of the mesh
        elFaces %A cell array containing electrode nodes
        colormap%colormap to use in plots
        scales  %the scales to multiply the estimates to get the plotted values
        hFig    %figure handle
        hUQFig  %Uncertainty qunatification figure handle
        hLineFig %Lineplot figure handle
        plotEl  %A flag: do we want to plot the electrodes? (This should be improved in future releases)
        title   %A title for the plot
        UQTitle %A title for the UQ plot
        linePlotTitle %A title for the lineplot
    end
    
    methods
        
        function obj = Plotter(g, H, elFaces)
            %Class constructor.
            %Input: g and H define the mesh we want to plot on. elFaces is
            %a cell array containing the nodes indices (rows of g) that
            %define electrodes.
            obj.g = g;
            obj.H = H;
            if nargin > 3
                obj.elFaces = elFaces;
                obj.plotEl = 1;
            else
                obj.plotEl = 0;
            end

            %Default values
            obj.colormap = jet;
            obj.scales = 1;
            obj.title = 'Conductivity';
            obj.UQTitle = 'Credible interval width';
            obj.linePlotTitle = 'Conductivity';
        end
        
        function plot2D(self, est, figh)
            %The basic plot function to plot a function in the imaging domain. Most likely
            %called by self.plot(...).
            
            %First check if we try to plot the electrodes without having
            %given the necessary variables:
            if self.plotEl && (isempty(self.elFaces))
                warning('Cannot plot electrodes since elfaces is not set');
                self.plotEl = 0;
            end
            
            %Set the figure we want to plot on:
            set(0, 'CurrentFigure', figh);
            clf;

            minval = min(est);
            maxval = max(est);

            if size(self.g,2) == 2
                fh = trisurf(self.H, self.g(:,1), self.g(:,2), est);%This is the line that does the main plotting
                view(2);%view from up so that it looks like 2D plot
                set(fh, 'edgecolor', 'none');%no edges
                set(fh, 'FaceColor', 'interp');%smooth colors between elements
                colormap(get(fh,'parent'), self.colormap);%set colormap
            elseif size(self.g,2) == 3
                tr = Otsu(est, 256, 0);%This is still kind of experimental, hopefully improved soon.
                sel = ~any(est(self.H)<tr,2);
                fh = tetramesh(self.H(sel,:),self.g);
            end
            if maxval > minval
                set(get(fh,'parent'),'CLim',[minval maxval]);
            else%Prevent error when plotting a constant sigma
                set(get(fh,'parent'),'CLim',[minval-1e-6 maxval+1e-6]);
            end
            axis('square')
            colorbar;%show the colorbar on the side of the plot
            %Plot the electrodes:
            if self.plotEl%We want to plot the electrodes as well
                hold on;
                for iel = 1:length(self.elFaces)
                    for it = 1:length(self.elFaces{iel})
                        plot(self.g(self.elFaces{iel}(it,:), 1), self.g(self.elFaces{iel}(it,:), 2), 'k-', 'LineWidth', 5);
                    end
                end
                hold off;
            end
            
        end

        function plotUQ(self, est)
            %Plot on the self.hUQFig figure. est should be a vector
            %defining the plotted function values on the nodes defined by g
            
            if isempty(self.hUQFig)
                self.hUQFig = figure();
            end
            est = est.*self.scales;
            self.plot2D(est, self.hUQFig);
            title(self.UQTitle);
                
        end

        function plot(self, est)
            %Plot on the self.hFig figure. est should be a vector
            %defining the plotted function values on the nodes defined by g

            if isempty(self.hFig)
                self.hFig = figure();
            end
            est = est.*self.scales;
            self.plot2D(est, self.hFig);
            title(self.title);
                
        end

        function lineplot(self, est, UQ, trueVal, p, t)
            %This function plots the estimate along a line through the
            %imaging domain. In addition, the credible interval (UQ) is
            %plotted. Optinally, a true value (trueVal) may be added as well for
            %comparison.
            %
            %Note: at the moment, this works only on 2D distributions!
            %
            %the plotting line is defined so that it starts from point p
            %and extends through the vector t. Default arguments for these
            %make the line start from the point with least x-value and at
            %y-value 0, and extend to the maximum x-value and y-value 0.

            if nargin < 5 || isempty(p)
                p = [min(self.g(:,1)); 0];
            end
            if nargin < 6 || isempty(t)
                t = [max(self.g(:,1))-min(self.g(:,1)); 0];
            end

            if isempty(self.hLineFig)
                self.hLineFig = figure();
            end
            est = est.*self.scales;
            UQ = UQ.*self.scales;
            trueVal = trueVal.*self.scales;

            points = p + [linspace(0,1,100)*t(1); linspace(0,1,100)*t(2)];
            points = points';
            PM = ForwardMesh1st.interpolatematrix2d(self.H, self.g, points);
            
            set(0, 'CurrentFigure', self.hLineFig);
            clf;
            xvals = linspace(min(self.g(:,1)), max(self.g(:,1)), 100);
            if isempty(trueVal)
                plot(xvals, PM*est, 'b-', xvals, PM*(est-UQ), 'r:', xvals, PM*(est+UQ), 'r:', 'LineWidth', 2);
                legend('Estimate', 'credible interval', '', 'Location', 'best');
            else
                plot(xvals, PM*est, 'b-', xvals, PM*(est-UQ), 'r:', xvals, PM*(est+UQ), 'r:', xvals, PM*trueVal, 'k-', 'LineWidth', 2);
                legend('Estimate', 'credible interval', '', 'True value', 'Location', 'best');
            end

            xlabel('x (m)');
            ylabel('\sigma (S/m)');
            set(gca, 'FontSize', 18);
            title(self.linePlotTitle);

        end

        
    end    
    
    
end