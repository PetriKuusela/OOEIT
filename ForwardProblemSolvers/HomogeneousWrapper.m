classdef HomogeneousWrapper < handle
%This is a wrapper class for EITFEM class used for computing the
%homogeneous estimates. It changes the input from being full sigma to just
%one value. Enables simultaneous estimation of contact impedances and the
%homogeneous estimate.
    
    properties
        solver  %This is the EITFEM object to be wrapped
        estInd  %in case of multiple estimates, which index to homogenize
        ng      %The number of nodes in the inverse mesh
    end
    
    methods
        
        function obj = HomogeneousWrapper(solver)
            %Class constructor. Input argument is the solver to be wrapped.
            obj.solver = solver;
            obj.ng = solver.fmesh.ng;
            obj.estInd = 1;
        end
        
        function vec = SolveForwardVec(self, hest)
            if isa(hest, 'EstimateVec')
                hest.estimates{self.estInd} = hest.estimates{self.estInd}*ones(self.ng,1);
                vec = self.solver.SolveForwardVec(hest);
            else
                vec = self.solver.SolveForwardVec(hest*ones(self.ng,1));
            end
        end
        
        function res = OptimizationFunction(self,hest)
            if isa(hest, 'EstimateVec')
                hest.estimates{self.estInd} = hest.estimates{self.estInd}*ones(self.ng,1);
                res = self.solver.OptimizationFunction(hest);
            else
                res = self.solver.OptimizationFunction(hest*ones(self.ng,1));
            end
        end
        
        function [Hess, grad] = GetHessAndGrad(self,hest)
            if isa(hest, 'EstimateVec')
                hest.estimates{self.estInd} = hest.estimates{self.estInd}*ones(self.ng,1);
                [Hess, grad] = self.solver.GetHessAndGrad(hest);
                grad.estimates{self.estInd} = sum(grad.estimates{self.estInd});
                Hess.estimates{self.estInd,self.estInd} = sum(sum(Hess.estimates{self.estInd, self.estInd}));
                for ii = 1:length(grad.estimates)
                    if ii < self.estInd
                        Hess.estimates{ii, self.estInd} = sum(Hess.estimates{ii,self.estInd}, 2);
                    elseif ii > self.estInd
                        Hess.estimates{self.estInd, ii} = sum(Hess.estimates{self.estInd,ii}, 1);
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
            if isa(hest, 'EstimateVec')
                hest.estimates{self.estInd} = hest.estimates{self.estInd}*ones(self.ng,1);
                self.solver.Plot(hest);
            else
                self.solver.Plot(hest*ones(self.ng,1));
            end
        end
        
        
    end
    
end