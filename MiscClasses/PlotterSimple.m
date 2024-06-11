classdef PlotterSimple < handle
%A Plotter class that can be assigned to inverse problem solver (e.g.
%SolverGN.m) to plot the estimate in between the iterations.
%
%This plotter plots the estimate as a simple lineplot. Useful e.g. for
%contact impedances.
    
    properties
        hFig     %figure handle
        hUQfig   %Uncertainty qunatification figure handle
        title   %A title for the basic figure
        UQTitle %A title for the UQ figure
    end
    
    methods

        function obj = PlotterSimple()
            %Default titles:
            obj.title = 'Contact impedances';
            obj.UQTitle = 'Contact impedance UQ';
        end

        function plot(self, est)
            %This function plots est as a simple lineplot into figure
            %self.hFig
            
            if isempty(self.hFig)
                self.hFig = figure();
            end
            set(0, 'CurrentFigure', self.hFig);
            plot(est);
            title(self.title);
                
        end

        function plotUQ(self, est)
            %This function plots est as a simple lineplot into figure
            %self.hUQFig
            
            if isempty(self.hUQfig)
                self.hUQfig = figure();
            end
            set(0, 'CurrentFigure', self.hUQfig);
            plot(est);
            title(self.UQTitle);
                
        end

        function lineplot(self, ~, ~, ~, ~, ~)
            %This function may be called by PlotterHolder class, but we do
            %not actually want to lineplot anything with this class.
        end
        
    end    
    
    
end