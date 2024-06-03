classdef PlotterSimple < handle
%A Plotter class that can be assigned to inverse problem solver (e.g.
%SolverLinesearch.m) to plot the estimate in between the iterations.
    
    properties
        fig     %figure handle
        UQfig   %Uncertainty qunatification figure handle
        title   %A title for the basic figure
        UQtitle %A title for the UQ figure
    end
    
    methods

        function obj = PlotterSimple()
            obj.title = 'Contact impedances';
            obj.UQtitle = 'Contact impedance UQ';
        end

        function plot(self, est)
            
            if isempty(self.fig)
                self.fig = figure();
            end
            set(0, 'CurrentFigure', self.fig);
            plot(est);
            title(self.title);
                
        end

        function plotUQ(self, est)
            
            if isempty(self.UQfig)
                self.UQfig = figure();
            end
            set(0, 'CurrentFigure', self.UQfig);
            plot(est);
            title(self.UQtitle);
                
        end

        
    end    
    
    
end