classdef PlotterHolder < handle
%A Plotter class that can be assigned to inverse problem solver (e.g.
%SolverLinesearch.m) to plot the estimate in between the iterations.
    
    properties
        plotters
        g       %The array of node co-ordinates
        H       %The elements of the mesh
        elfaces %A cell array containing electrode nodes
        cm      %colormap
        scales  %the scales to multiply the estimates to get the plotted values
        fig     %figure handle(s)
        plottypes
        UQfig   %Uncertainty qunatification figure handle(s)
        linefig %Lineplot figure handle(s)
        plotel  %A flag: do we want to plot the electrodes?
    end
    
    methods
        
        function obj = PlotterHolder(plotters)
            %Class constructor.
            obj.plotters = plotters;
        end
        
        function plot(self, est)
            
            
            if ~isa(est, 'Estimate_vec')
                self.plotters{1}.plot(est);
            else
                for ii = 1:length(est.estimates)
                    if length(self.plotters) >= ii && ~isempty(self.plotters{ii})
                        self.plotters{ii}.plot(est.estimates{ii});
                    end
                end
            end
                
        end
        
        function plotUQ(self, est)
            
            if ~isa(est, 'Estimate_vec')
                self.plotters{1}.plotUQ(est);
            else
                for ii = 1:length(est.estimates)
                    if length(self.plotters) >= ii && ~isempty(self.plotters{ii})
                        self.plotters{ii}.plotUQ(est.estimates{ii});
                    end
                end
            end
        end

        function lineplot(self, est, UQ, trueval, p, t)

            if nargin < 5 || isempty(p)
                p = [min(self.g(:,1)); 0];
            end
            if nargin < 6 || isempty(t)
                t = [max(self.g(:,1))-min(self.g(:,1)); 0];
            end

            if ~isa(est, 'Estimate_vec')
                self.plotters{1}.lineplot(est, UQ, trueval, p, t);
            else
                for ii = 1:length(est.estimates)
                    if length(self.plotters) >= ii && ~isempty(self.plotters{ii})
                        self.plotters{ii}.lineplot(est.estimates{ii}, UQ.estimates{ii}, trueval.estimates{ii}, p{ii}, t{ii});
                    end
                end
            end

        end
        
        
        
    end    
    
    
end