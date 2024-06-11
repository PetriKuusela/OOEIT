classdef EITFEMComplex < EITFEM
%This class computes complex CEM solutions and takes the complex
%conductivity as an Estimate_vec object with real and imag parts in their
%own cells.
    properties

        epsilonInd
        omega

    end


    methods

        function obj = EITFEMComplex(fmesh)
            obj = obj@EITFEM(fmesh);
            obj.epsilonInd = 2;
            obj.omega = 1;
        end

        function elval = SolveForward(self, est, omega)
            if nargin < 3 || isempty(omega)
                omega = self.omega;
            end
            sigma = est.estimates{self.sigmaInd} + 1i*omega*est.estimates{self.epsilonInd};
            elval = SolveForward@EITFEM(self, sigma);
        end

        function vec = SolveForwardVec(self, est, omega)
            if nargin < 3 || isempty(omega)
                omega = self.omega;
            end
            vec = self.SolveForward(est, omega);
            vec = vec(:);
            if ~isempty(self.mIncl)
                vec = vec(self.mIncl);
            end
            vec = [real(vec); imag(vec)];
            vec = vec + self.eps;
            
        end

        function SetInvGamma(self, constNoise, relNoise)
            if nargin < 3 || isempty(relNoise)
                relNoise = [0; 0];
            end
            
            if length(constNoise) == 1
                constNoise = [constNoise; constNoise];
            end
            if length(relNoise) == 1
                relNoise = [relNoise; relNoise];
            end

            %Are the noise levels determined based on I or U?
            if strcmp(self.mode, 'current')
                meas = self.Uel;
            elseif strcmp(self.mode, 'potential')
                meas = self.Iel;
            else
                error(['Unrecognized solver mode: ' self.mode]);
            end
            nmeas = length(meas);
            
            var_meas = zeros(nmeas,1);
            %Calculating the model variance of the noise
            select = 1:nmeas/2;
            var_meas(select) = (constNoise(1)*(max(abs(meas(select)))-min(abs(meas(select)))))^2;
            var_meas(select) = var_meas(select) + (relNoise(1)*(abs(meas(select)))).^2;
            %imag:
            select = nmeas/2+1:nmeas;
            var_meas(select) = (constNoise(2)*(max(abs(meas(select)))-min(abs(meas(select)))))^2;
            var_meas(select) = var_meas(select) + (relNoise(2)*(abs(meas(select)))).^2; %Assume no cross-correlations, hence a diagonal covariance mat:
            Gamma = diag(var_meas(:));
            self.InvGamma = sparse(inv(Gamma));
            self.Ln = chol(self.InvGamma);%Store both invGamma and it's Cholesky
        end

        function [Hess, grad] = GetHessAndGrad(self, est, omega)
            if nargin < 3 || isempty(omega)
                omega = self.omega;
            end

            HC = cell(length(est.estimates));%Initialize the Hess and gradient cell arrays for the return objects
            gC = cell(length(est.estimates),1);

            sigma = est.estimates{self.sigmaInd} + 1i*omega*est.estimates{self.epsilonInd};
            fRes = self.SolveForwardVec(est, omega);
            if strcmp(self.mode, 'potential')
                diff = (fRes - self.Iel);
            elseif strcmp(self.mode, 'current')
                diff = (fRes - self.Uel);
            else
                error(['Unrecognized solver mode: ' self.mode]);
            end
            J = self.Jacobian(sigma, 1);
            J = [real(J) -omega*imag(J); imag(J) omega*real(J)];
            if length(self.scales) == 2
                J(:,1:end/2) = J(:,1:end/2)*self.scales(1);
                J(:,end/2+1:end) = J(:,end/2+1:end)*self.scales(2);
            else
                J = J*diag(self.scales);%The scaling done before solving the FEM has to be taken into account here.
            end
            Hess = J'*self.InvGamma*J;%First order approximation of the forward problem
            grad = J'*self.InvGamma*diff;
            HC{self.sigmaInd, self.sigmaInd} = Hess(1:end/2,1:end/2);
            HC{self.epsilonInd, self.epsilonInd} = Hess(end/2+1:end,end/2+1:end);
            HC{self.sigmaInd, self.epsilonInd} = Hess(1:end/2,end/2+1:end);
            gC{self.sigmaInd} = grad(1:end/2);
            gC{self.epsilonInd} = grad(end/2+1:end);
            Hess = EstimateHess(HC);
            grad = EstimateVec(gC);

        end

        function sigma = PreProcessSigma(self, sigma)
            
            if length(self.sigmaMin) == 1
                sigMin = [self.sigmaMin; self.sigmaMin];
            else
                sigMin = self.sigmaMin;
            end
            if length(self.scales) == 2
                rsigma = self.scales(1)*real(sigma);
                isigma = self.scales(2)*imag(sigma);
                rsigma(rsigma<sigMin(1)) = sigMin(1);
                isigma(isigma<sigMin(2)) = sigMin(2);
                sigma = rsigma + 1i*isigma;
            elseif length(self.scales) > 1
                rsigma = self.scales(1:end/2).*real(sigma);
                isigma = self.scales(end/2+1:end).*imag(sigma);
                rsigma(rsigma<sigMin(1)) = sigMin(1);
                isigma(isigma<sigMin(2)) = sigMin(2);
                sigma = rsigma + 1i*isigma;
            else
                sigma = self.scales*sigma;
                sigma(real(sigma)<sigMin(1)) = sigMin(1) + 1i*imag(sigma(real(sigma)<sigMin(1)));
                sigma(imag(sigma)<sigMin(2)) = real(sigma(imag(sigma)<sigMin(2))) + 1i*sigMin(2);
            end
        
        end

    end

end