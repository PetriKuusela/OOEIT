classdef EITFEMNLD < EITFEM
%This class is used for computing non-linear difference reconstructions.
%This has to be used in conjunction with normal EITFEM object, which solves
%CEM for the initial distribution, and this object solves the CEM for the
%final distribution = initial distribution + change.

    properties
        sigmaInitialInd
    end

    methods

        function obj = EITFEMNLD(fmesh)
            obj = obj@EITFEM(fmesh);
            obj.sigmaInitialInd = 1;
            obj.sigmaInd = 2;
        end

        function elval = SolveForwardVec(self, est)
        
            elval = SolveForwardVec@EITFEM(self, est.estimates{self.sigmaInd} + est.estimates{self.sigmaInitialInd});
            
        end

        function [Hess, grad] = GetHessAndGrad(self, est)

            [Hess, grad] = GetHessAndGrad@EITFEM(self, est);
            Hess.estimates{self.sigmaInitialInd, self.sigmaInd} = Hess.estimates{self.sigmaInd, self.sigmaInd};
            Hess.estimates{self.sigmaInitialInd, self.sigmaInitialInd} = Hess.estimates{self.sigmaInd, self.sigmaInd};
            grad.estimates{self.sigmaInitialInd} = grad.estimates{self.sigmaInd};

        end


    end


end