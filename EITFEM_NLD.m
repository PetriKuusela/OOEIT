classdef EITFEM_NLD < EITFEM
%This class is used for computing non-linear difference reconstructions.
%This has to be used in conjunction with normal EITFEM object, which solves
%CEM for the initial distribution, and this object solves the CEM for the
%final distribution = initial distribution + change.

    properties
        sigma_iind
    end

    methods

        function obj = EITFEM_NLD(fmesh)
            obj = obj@EITFEM(fmesh);
            obj.sigma_iind = 1;
            obj.sigmaind = 2;
        end

        function elval = SolveForwardVec(self, est)
        
            elval = SolveForwardVec@EITFEM(self, est.estimates{self.sigmaind} + est.estimates{self.sigma_iind});
            
        end

        function [Hess, grad] = GetHessAndGrad(self, est)

            [Hess, grad] = GetHessAndGrad@EITFEM(self, est);
            Hess.estimates{self.sigma_iind, self.sigmaind} = Hess.estimates{self.sigmaind, self.sigmaind};
            Hess.estimates{self.sigma_iind, self.sigma_iind} = Hess.estimates{self.sigmaind, self.sigmaind};
            grad.estimates{self.sigma_iind} = grad.estimates{self.sigmaind};

        end



    end


end