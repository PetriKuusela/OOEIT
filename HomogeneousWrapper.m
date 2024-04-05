classdef HomogeneousWrapper < handle
%This is a wrapper class for EITFEM class used for computing the
%homogeneous estimates. It changes the input from being full sigma to just
%one value. Enables simultaneous estimation of contact impedances and the
%homogeneous estimate.
    
    properties
        solver  %This is the EITFEM object to be wrapped
        estind  %in case of multiple estimates, which index to homogenize
        ng      %The number of nodes in the inverse mesh
    end
    
    methods
        
        function obj = HomogeneousWrapper(solver)
            %Class constructor. Input argument is the solver to be wrapped.
            obj.solver = solver;
            obj.ng = solver.fmesh.ng;
            obj.estind = 1;
        end
        function vec = SolveForwardVec(self, hest)
            if isa(hest, 'Estimate_vec')
                hest.estimates{self.estind} = hest.estimates{self.estind}*ones(self.ng,1);
                vec = self.solver.SolveForwardVec(hest);
            else
                vec = self.solver.SolveForwardVec(hest*ones(self.ng,1));
            end
        end
        function res = OptimizationFunction(self,hest)
            if isa(hest, 'Estimate_vec')
                hest.estimates{self.estind} = hest.estimates{self.estind}*ones(self.ng,1);
                res = self.solver.OptimizationFunction(hest);
            else
                res = self.solver.OptimizationFunction(hest*ones(self.ng,1));
            end
        end
        function [Hess, grad] = GetHessAndGrad(self,hest)
            if isa(hest, 'Estimate_vec')
                hest.estimates{self.estind} = hest.estimates{self.estind}*ones(self.ng,1);
                [Hess, grad] = self.solver.GetHessAndGrad(hest);
                grad.estimates{self.estind} = sum(grad.estimates{self.estind});
                Hess.estimates{self.estind,self.estind} = sum(sum(Hess.estimates{self.estind, self.estind}));
                for ii = 1:length(grad.estimates)
                    if ii < self.estind
                        Hess.estimates{ii, self.estind} = sum(Hess.estimates{ii,self.estind}, 2);
                    elseif ii > self.estind
                        Hess.estimates{self.estind, ii} = sum(Hess.estimates{self.estind,ii}, 1);
                    end
                end
            else
                [Hess, grad] = self.solver.GetHessAndGrad(hest*ones(self.ng,1));
                grad = sum(grad);
                Hess = sum(sum(Hess(1:end/2,1:end/2)));
            end
        end
        function Plot(self,hest)
            %This can be used to plot the fitting of the data during the
            %solving of the homogeneous estimate.
            if isa(hest, 'Estimate_vec')
                hest.estimates{self.estind} = hest.estimates{self.estind}*ones(self.ng,1);
                self.solver.Plot(hest);
            else
                self.solver.Plot(hest*ones(self.ng,1));
            end
        end
        
        
    end
    
end