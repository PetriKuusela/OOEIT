classdef PriorSmoothness < handle
%PriorSmoothness class contains a smoothness prior and the required functions
%for it to be called by an inverse problem solver -class (e.g.
%SolverLinesearch.m).
%
%Author: Petri Kuusela 5.4.2024
    
    
    properties
        invCov  %The inverse of the covariance matrix
        mean    %The expectation value of the estimate (can be a single value or a vector with same length as ginv)
        L       %chol(invcov), used only to draw samples from the prior
        regularization_c %A regularization constant to make sure the covariance matrix is not poorly conditioned
        common_c %A constant added to all the elements of the covariance matrix. Larger values allow greater variation in the mean value of the estimate
    end
    
    
    methods
        
        function obj = PriorSmoothness(ginv, corlen, var, mean)
            %Class constructor.
            %Input: ginv = the nodes of the reconstruction mesh
            %       corlen = the correlation length (length where the
            %               cross-covariance has decreased to 1%)
            %       var = the variance of the values of the estimate
            %       mean = the expectation value of the estimate.
            
            obj.regularization_c = 1e-4;
            obj.common_c = 1e-1;
            obj.SetCovMat(ginv, corlen, var);
            obj.mean = mean;
            obj.L = chol(obj.invCov);
        end
        
        function SetCovMat(self, g, corlen, var)
            %Compute the inverse covariance matrix given the input
            %parameters. 
            %Input: g = node coordinates
            %       corlen = spatial correlation length (length where the
            %               cross-covariance has decreased to 1%)
            %       var = the variance of the estimate values
            
            ng = size(g,1);
            b = corlen./sqrt(2*log(100));
            c = self.regularization_c*var;
            a = var-c;
            xmat = repmat(g(:,1),1,ng);
            ymat = repmat(g(:,2),1,ng);
            zmat = 0;
            if size(g,2) > 2%if the mesh is 3D
                zmat = repmat(g(:,2),1,ng);
            end
            cov = a*exp(-0.5*((xmat-xmat').^2+(ymat-ymat').^2+(zmat-zmat').^2)/b^2);
            cov = cov + diag(c(1)*ones(size(cov,1),1)) + self.common_c*var(1)*ones(size(cov,1));
            self.invCov = inv(cov);
        end
        
        function res = OptimizationFunction(self, sigest)
            %This function is called by the inverse problem solver, and it
            %gives the function value to be minimized.
            res = 0.5*(sigest-self.mean)'*self.invCov*(sigest-self.mean);
        end
        
        function [Hess, grad] = GetHessAndGrad(self, sigest)
            %This function is called by the inverse problem solver, and it
            %gives the Hess-matrix and gradient of the optimization
            %function.
            grad = self.invCov*(sigest-self.mean);
            Hess = self.invCov;
        end

        function samples = DrawSamples(self, n)
            vec = randn(size(self.L,1),n);
            if length(self.mean) == 1
                samples = self.L\vec + self.mean*ones(size(self.L,1),n);
            else
                samples = self.L\vec + self.mean*ones(1,n);
            end

        end

        
        
    end
    
    
    
end