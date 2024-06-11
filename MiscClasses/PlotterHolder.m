classdef PlotterHolder < handle
%A Plotter class that can be assigned to inverse problem solver (e.g.
%SolverGN.m) to plot the estimate in between the iterations.
%
%The class contains a cell array of plotters which do the plotting. Cell
%array indexing corresponds to the indexing on the estimate variable i.e.
%the cell array indexing on EstimateVec.estimates. For estimates not to be
%plotted, empty elements may be left into the cell array PlotterHolder.plotters.
    
    properties
        plotters%This contains a cell array of plotters that do the plotting
        %g       %The array of node co-ordinates
        %H       %The elements of the mesh
        %elfaces %A cell array containing electrode nodes
        %cm      %colormap
        %scales  %the scales to multiply the estimates to get the plotted values
        %fig     %figure handle(s)
        %plottypes
        %UQfig   %Uncertainty qunatification figure handle(s)
        %linefig %Lineplot figure handle(s)
        %plotel  %A flag: do we want to plot the electrodes?
    end
    
    methods
        
        function obj = PlotterHolder(plotters)
            %Class constructor.
            obj.plotters = plotters;
        end
        
        function plot(self, est)
            
            
            if ~isa(est, 'EstimateVec')
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
            
            if ~isa(est, 'EstimateVec')
                self.plotters{1}.plotUQ(est);
            else
                for ii = 1:length(est.estimates)
                    if length(self.plotters) >= ii && ~isempty(self.plotters{ii})
                        self.plotters{ii}.plotUQ(est.estimates{ii});
                    end
                end
            end
        end

        function lineplot(self, est, UQ, trueVal, p, t)

            if ~isa(est, 'EstimateVec')
                self.plotters{1}.lineplot(est, UQ, trueVal, p, t);
            else
                for ii = 1:length(est.estimates)
                    if length(self.plotters) >= ii && ~isempty(self.plotters{ii})
                        self.plotters{ii}.lineplot(est.estimates{ii}, UQ{ii}, trueVal{ii}, p{ii}, t{ii});
                    end
                end
            end

        end
        
        
        
    end    
    
    
end