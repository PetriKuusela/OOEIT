classdef Plotter3D < handle
%A Plotter class that can be assigned to inverse problem solver (e.g.
%SolverGN.m) to plot the estimate in between the iterations.
    
    properties
        g       %The array of node co-ordinates
        H       %The elements of the mesh
        edges   %node pairs describing edgest of elements
        elFaces %A cell array containing electrode nodes
        colormap%colormap to use in plots
        scales  %the scales to multiply the estimates to get the plotted values
        hFig    %figure handle
        hUQFig  %Uncertainty quantification figure handle
        hLineFig %Lineplot figure handle
        plotEl  %A flag: do we want to plot the electrodes? (This should be improved in future releases)
        title   %A title for the plot
        UQTitle %A title for the UQ plot
        linePlotTitle %A title for the lineplot
        basex   %middle coordinates for intersectorplanes to use in plot3D
        basey
        basez
    end
    
    methods
        
        function obj = Plotter3D(g, H, elFaces)
            %Class constructor.
            %Input: g and H define the mesh we want to plot on. elFaces is
            %a cell array containing the nodes indices (rows of g) that
            %define electrodes.
            obj.g = g;
            obj.H = H;
            obj.basex = 0.5*(min(g(:,1)) + max(g(:,1)));
            obj.basey = 0.5*(min(g(:,2)) + max(g(:,2)));
            obj.basez = 0.5*(min(g(:,3)) + max(g(:,3))) +1e-5; %fix this
            obj.edges = obj.FindEdges(H);
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

        function plot3D(self, est, figh, x, y, z)

            
            %Set the figure we want to plot on:
            set(0, 'CurrentFigure', figh);
            clf;


            hold on
            for ii = 1:length(x)
                self.PlotSurface(est, x(ii), 'x');
            end
            for ii = 1:length(y)
                self.PlotSurface(est, y(ii), 'y');
            end
            for ii = 1:length(z)
                self.PlotSurface(est, z(ii), 'z');
            end

            
            axis equal;
            colorbar;
        end

        function plotUQ(self, est)
            %Plot on the self.hUQFig figure. est should be a vector
            %defining the plotted function values on the nodes defined by g
            
            if isempty(self.hUQFig)
                self.hUQFig = figure();
            end
            est = est.*self.scales;
            self.plot3D(est, self.hUQFig, self.basex, self.basey, self.basez);
            title(self.UQTitle);
                
        end

        function plot(self, est)
            %Plot on the self.hFig figure. est should be a vector
            %defining the plotted function values on the nodes defined by g

            if isempty(self.hFig)
                self.hFig = figure();
            end
            est = est.*self.scales;
            self.plot3D(est, self.hFig, self.basex, self.basey, self.basez);
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

        function ed = FindEdges(self, H)
            ee = 6;
            ed = zeros(ee*size(H,1),2);
            for ii = 1:size(H,1)
                ed(ee*(ii-1)+1,:) = H(ii,[1 2]);
                ed(ee*(ii-1)+2,:) = H(ii,[1 3]);
                ed(ee*(ii-1)+3,:) = H(ii,[1 4]);
                ed(ee*(ii-1)+4,:) = H(ii,[2 3]);
                ed(ee*(ii-1)+5,:) = H(ii,[2 4]);
                ed(ee*(ii-1)+6,:) = H(ii,[3 4]);
            end
            ed = [ed(:,1) ed(:,2); ed(:,2) ed(:,1)];
            del = ed(:,2) < ed(:,1);
            ed(del,:) = [];
            ed = unique(ed, 'rows');
        end

        function PlotSurface(self, vals, c, type)
            if strcmp(type, 'x')
                t = [self.g(self.edges(:,1),1)-c self.g(self.edges(:,2),1)-c];
            elseif strcmp(type, 'y')
                t = [self.g(self.edges(:,1),2)-c self.g(self.edges(:,2),2)-c];
            elseif strcmp(type, 'z')
                t = [self.g(self.edges(:,1),3)-c self.g(self.edges(:,2),3)-c];
            else
                error(['Unknown type: ' type]);
            end
            sel = t(:,1).*t(:,2)<0; %select elements where  t(ii,1) and t(ii,2) have opposite signs
            
            ttot = sum(abs(t(sel, :)),2);
            evals = abs(t(sel,1)).*vals(self.edges(sel,2)) + abs(t(sel,2)).*vals(self.edges(sel,1));
            evals = evals./ttot;

            gsurfx = (abs(t(sel,1)).*self.g(self.edges(sel,2),1) + abs(t(sel,2)).*self.g(self.edges(sel,1),1))./ttot;
            gsurfy = (abs(t(sel,1)).*self.g(self.edges(sel,2),2) + abs(t(sel,2)).*self.g(self.edges(sel,1),2))./ttot;
            gsurfz = (abs(t(sel,1)).*self.g(self.edges(sel,2),3) + abs(t(sel,2)).*self.g(self.edges(sel,1),3))./ttot;

            gsurf = [gsurfx gsurfy gsurfz];

            if strcmp(type, 'x')
                Hsurf = delaunay(gsurfy, gsurfz);
            elseif strcmp(type, 'y')
                Hsurf = delaunay(gsurfx, gsurfz);
            elseif strcmp(type, 'z')
                Hsurf = delaunay(gsurfx, gsurfy);
            else
                error(['Unknown type: ' type]);
            end
            h = trisurf(Hsurf, gsurf(:,1), gsurf(:,2), gsurf(:,3), evals);
            set(h, 'edgecolor', 'none');

        end

        
    end    
    
    
end