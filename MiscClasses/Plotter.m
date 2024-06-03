classdef Plotter < handle
%A Plotter class that can be assigned to inverse problem solver (e.g.
%SolverLinesearch.m) to plot the estimate in between the iterations.
    
    properties
        g       %The array of node co-ordinates
        H       %The elements of the mesh
        elfaces %A cell array containing electrode nodes
        cm      %colormap
        scales  %the scales to multiply the estimates to get the plotted values
        fig     %figure handle(s)
        UQfig   %Uncertainty qunatification figure handle(s)
        linefig %Lineplot figure handle(s)
        plotel  %A flag: do we want to plot the electrodes?
        title   %A title for the plot
        UQtitle %A title for the UQ plot
        lineplottitle %A title for the lineplot
    end
    
    methods
        
        function obj = Plotter(g, H, elfaces)
            %Class constructor.
            %Input: g and H define the mesh we want to plot on
            obj.g = g;
            obj.H = H;
            obj.cm = jet;
            if nargin > 3
                obj.elfaces = elfaces;
            end
            obj.plotel = 0;
            obj.scales = 1;
            obj.title = 'Conductivity';
            obj.UQtitle = 'Credible interval width';
            obj.lineplottitle = 'Conductivity';
        end
        
        function plot2D(self, est, figh)
            %The basic plot function to plot sigma.
            
            %First check if we try to plot the electrodes without having
            %given the necessary variables:
            if self.plotel && (isempty(self.elfaces))
                warning('Cannot plot electrodes since elfaces is not set');
                self.plotel = 0;
            end
            
            %Set the figure we want to plot on:
            set(0, 'CurrentFigure', figh);
            clf;

            minval = min(est);
            maxval = max(est);

            if size(self.g,2) == 2
                fh = trisurf(self.H, self.g(:,1), self.g(:,2), est);%This is the line that does the main plotting
                view(2);
                set(fh, 'edgecolor', 'none');%no edges
                set(fh, 'FaceColor', 'interp');%smooth colors between elements
                colormap(get(fh,'parent'), self.cm);%set colormap
            elseif size(self.g,2) == 3
                tr = Otsu(est, 256, 0);
                sel = ~any(est(self.H)<tr,2);
                fh = tetramesh(self.H(sel,:),self.g);
            end
            if maxval > minval
                set(get(fh,'parent'),'CLim',[minval maxval]);
            else%Prevent error when plotting a constant sigma
                set(get(fh,'parent'),'CLim',[minval-1e-6 maxval+1e-6]);
            end
            axis('square')
            colorbar;
            %Plot the electrodes:
            if self.plotel
                hold on;
                for iel = 1:length(self.elfaces)
                    for it = 1:length(self.elfaces{iel})
                        plot(self.g(self.elfaces{iel}(it,:), 1), self.g(self.elfaces{iel}(it,:), 2), 'k-', 'LineWidth', 5);
                    end
                end
                hold off;
            end
            
        end

        function plotUQ(self, est)
            
            if isempty(self.UQfig)
                self.UQfig = figure();
            end
            est = est.*self.scales;
            self.plot2D(est, self.UQfig);
            title(self.UQtitle);
                
        end

        function plot(self, est)
            
            if isempty(self.fig)
                self.fig = figure();
            end
            est = est.*self.scales;
            self.plot2D(est, self.fig);
            title(self.title);
                
        end

        function lineplot(self, est, UQ, trueval, p, t)

            if nargin < 5 || isempty(p)
                p = [min(self.g(:,1)); 0];
            end
            if nargin < 6 || isempty(t)
                t = [max(self.g(:,1))-min(self.g(:,1)); 0];
            end

            if isempty(self.linefig)
                self.linefig = figure();
            end
            est = est.*self.scales;
            UQ = UQ.*self.scales;
            trueval = trueval.*self.scales;

            points = p + [linspace(0,1,100)*t(1); linspace(0,1,100)*t(2)];
            points = points';
            PM = ForwardMesh1st.interpolatematrix2d(self.H, self.g, points);
            
            set(0, 'CurrentFigure', self.linefig);
            clf;
            xvals = linspace(min(self.g(:,1)), max(self.g(:,1)), 100);
            if isempty(trueval)
                plot(xvals, PM*est, 'b-', xvals, PM*(est-UQ), 'r:', xvals, PM*(est+UQ), 'r:', 'LineWidth', 2);
                legend('Estimate', 'credible interval', '', 'Location', 'best');
            else
                plot(xvals, PM*est, 'b-', xvals, PM*(est-UQ), 'r:', xvals, PM*(est+UQ), 'r:', xvals, PM*trueval, 'k-', 'LineWidth', 2);
                legend('Estimate', 'credible interval', '', 'True value', 'Location', 'best');
            end

            xlabel('x (m)');
            ylabel('\sigma (S/m)');
            set(gca, 'FontSize', 18);
            title(self.lineplottitle);

        end

        
    end    
    
    
end