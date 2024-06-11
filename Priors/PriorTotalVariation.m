classdef PriorTotalVariation < handle
%PriorTotalVariation class contains a TV prior and the required functions
%for it to be called by an inverse problem solver -class (e.g.
%SolverGN.m).
%
%Author: Petri Kuusela 16.11.2021
    
    properties
        gradphi     %The gradients of all basis functions. A cell array with each cell containing gradients inside one element
        H           %The elements of the reconstruction mesh
        g           %The nodes of the reconstruction mesh
        nH          %Number of elements
        ng          %Number of nodes
        alpha       %The "strength" coefficient of TV
        beta        %The smoothing parameter
        areas       %A vector containing the area of each element
        sigmaInd    %An array specifying which elements of EstimateVec should be taken into account
    end
    
    methods
        
        function obj = PriorTotalVariation(g, H, ec, ef, beta)
            %Class constructor.
            %Inputs: g = the nodes of the reconstruction mesh
            %        H = the elements of the reconstruction mesh
            %        ec = expected change of the estimate 
            %        ef(optional) = extra factor which multiplies alpha
            %        beta(optional) = the smoothing constant
            %NOTE: if you want to set alpha yourself, you can change it
            %directly after initiating the object.
            
            if nargin < 4 || isempty(ef)
                ef = 1;
            end
            if nargin < 5 || isempty(beta)
                beta = 1e-4;
            end
            
            obj.g = g;
            obj.H = H;
            mesh_a = (max(max(g(:,1)))-min(min(g(:,1))))*(max(max(g(:,2)))-min(min(g(:,2))))/size(H,1);%roughly the element area of the mesh
            mesh_c = sqrt(2*mesh_a);%roughly the element size of mesh
            egrad = ec./mesh_c;%expected gradient
            tvp = 0.975;
            obj.alpha = -ef.*log(1-tvp)./(mesh_a*egrad);
            obj.beta = beta;
            obj.sigmaInd = 1;
            
            obj.ng = size(g, 1);
            obj.nH = size(H, 1);
            [obj.gradphi, obj.areas] = obj.ComputeGrads2D(g,H);
            
        end
        
        function res = OptimizationFunction(self, est)
            %This is the function called by the inverse problem solver to
            %get the value of the function which is to be minimized.
            if isa(est, 'EstimateVec')
                res = 0;
                for ii = 1:length(self.sigmaInd)
                    gradsigma = self.ComputeGrads(est.estimates{self.sigmaInd(ii)});
                    res = res + self.alpha(ii)*sum(self.areas.*sqrt(sum(abs(gradsigma).^2, 2)+self.beta));
                end
            else
                gradsigma = self.ComputeGrads(est);
                res = self.alpha*sum(self.areas.*sqrt(sum(abs(gradsigma).^2, 2)+self.beta));
            end
        end
        
        function [Hess, grad] = GetHessAndGrad(self, est)
            %Compute the Hess-matrix and gradient of the optimization
            %function. This function is called by the inverse problem
            %solver.

            if isa(est, 'EstimateVec')
                HC = cell(length(est.estimates));
                gC = cell(length(est.estimates),1);
            end
            
            for io = 1:length(self.sigmaInd)
                grad = zeros(self.ng,1);%initialize arrays for gradient and Hess-matrix
                Hvals = zeros(9*self.nH,1);%Hess matrix is collected in sparse form
                HindsI = zeros(9*self.nH,1);
                HindsJ = zeros(9*self.nH,1);
                if isa(est, 'EstimateVec')
                    gradsigma = self.ComputeGrads(est.estimates{self.sigmaInd(io)});
                else
                    gradsigma = self.ComputeGrads(est);
                end

                II = 1;%index for collecting the Hess-values in sparce matrix
                %This block calculates the gradient of the TV functional
                for ii=1:self.nH
                    tdgrads = self.gradphi{ii};%The gradients of the basis functions in element ii
                    tinds = self.H(ii,:);%node indices that are the verices of element ii
                    tdsigma = gradsigma(ii,:);%The gradient of sigma in element ii
                    gnorm = 1/sqrt(sum(tdsigma.^2) + self.beta);%1/sqrt part of the gradient
                    tdot = tdgrads*tdsigma';%The dot product of gradients of basis functions and gradient of sigma
                    grad(tinds) = grad(tinds) + self.alpha(io)*gnorm.*tdot.*self.areas(ii);%Add the terms to the relevant gradient elements
                end

                %This block calculates the Hess-matrix of the TV functional
                for ii=1:self.nH
                    tdgrads = self.gradphi{ii};%The gradients of the basis functions in element ii
                    tinds = self.H(ii,:);%node indices that are the vertices of element ii
                    tdsigma = gradsigma(ii,:);%The gradient of sigma in element ii
                    gnorm = (sum(tdsigma.^2) + self.beta)^(-0.5);%These are parts of the Hessian
                    gnorm2 = (sum(tdsigma.^2) + self.beta)^(-1.5);
                    
                    for jj=1:3
                        for kk=1:3
                            phii = tdgrads(jj,:);%grad? of basis function_i
                            phij = tdgrads(kk,:);%grad? of basis function_j
                            f = phii*phij'*gnorm;
                            h = -(tdsigma*phii'*gnorm2*tdsigma*phij');
                            Hvals(II) = self.alpha(io)*(f + h).*self.areas(ii);
                            HindsI(II) = tinds(jj);
                            HindsJ(II) = tinds(kk);
                            II = II + 1;
                        end
                    end
                end

                if isa(est, 'EstimateVec')
                    HC{self.sigmaInd(io), self.sigmaInd(io)} = accumarray([HindsI HindsJ], Hvals);
                    gC{self.sigmaInd(io)} = grad;
                else
                    Hess = accumarray([HindsI HindsJ], Hvals);
                end

            end

            if isa(est,'EstimateVec')
                Hess = EstimateHess(HC);
                grad = EstimateVec(gC);
            end
            
        end
        
        function gradsigma = ComputeGrads(self, sigest)
            %Compute the gradients of the sigma, using the pre-computed
            %self.gradphi.
            N = size(self.gradphi,1);
            gradsigma = zeros(N,2);
            for ii=1:N
                sigmas = sigest(self.H(ii,:));
                grads = self.gradphi{ii};
                gradsigma(ii,:) = sigmas(1)*grads(1,:) +  sigmas(2)*grads(2,:) +  sigmas(3)*grads(3,:);
            end
        end
        
        function [gradphi, areas] = ComputeGrads2D(self, g,H)
            % Computes the spatial gradients of the linear basis functions
            % Also calculates the area of each element.
            N = size(H,1);
            gN = size(g,1);
            gradphi = cell(N,1);
            areas = zeros(N,1);
            L = [-1 1 0; -1 0 1];
            R = sparse(2*N,gN);
            for ii=1:N
                X = g(H(ii,:),:);
                Jt = L*X;
                areas(ii) = 1/2*abs(det(Jt));
                grads = Jt\L;
                grads = grads';
                gradphi{ii} = grads;
            end

        end

        
        
    end 
    
end