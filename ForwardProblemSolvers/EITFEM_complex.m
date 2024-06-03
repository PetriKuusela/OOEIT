classdef EITFEM_complex < EITFEM
%This class computes complex CEM solutions and takes the complex
%conductivity as an Estimate_vec object with real and imag parts in their
%own cells.
    properties

        epsilonind
        omega

    end


    methods

        function obj = EITFEM_complex(fmesh)
            obj = obj@EITFEM(fmesh);
            obj.epsilonind = 2;
            obj.omega = 1;
        end

        function elval = SolveForward(self, est, omega)
            if nargin < 3 || isempty(omega)
                omega = self.omega;
            end
            sigma = est.estimates{self.sigmaind} + 1i*omega*est.estimates{self.epsilonind};
            elval = SolveForward@EITFEM(self, sigma);
        end

        function vec = SolveForwardVec(self, est, omega)
            if nargin < 3 || isempty(omega)
                omega = self.omega;
            end
            vec = self.SolveForward(est, omega);
            vec = vec(:);
            if ~isempty(self.vincl)
                vec = vec(self.vincl);
            end
            vec = [real(vec); imag(vec)];
            vec = vec + self.eps;
            
        end

        function SetInvGamma(self, meas_noise_coef_e, meas_noise_coef2)
            if nargin < 3 || isempty(meas_noise_coef2)
                meas_noise_coef2 = [0; 0];
            end
            
            if length(meas_noise_coef_e) == 1
                meas_noise_coef_e = [meas_noise_coef_e; meas_noise_coef_e];
            end
            if length(meas_noise_coef2) == 1
                meas_noise_coef2 = [meas_noise_coef2; meas_noise_coef2];
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
            var_meas(select) = (meas_noise_coef_e(1)*(max(abs(meas(select)))-min(abs(meas(select)))))^2;
            var_meas(select) = var_meas(select) + (meas_noise_coef2(1)*(abs(meas(select)))).^2;
            %imag:
            select = nmeas/2+1:nmeas;
            var_meas(select) = (meas_noise_coef_e(2)*(max(abs(meas(select)))-min(abs(meas(select)))))^2;
            var_meas(select) = var_meas(select) + (meas_noise_coef2(2)*(abs(meas(select)))).^2; %Assume no cross-correlations, hence a diagonal covariance mat:
            Gamma_n = diag(var_meas(:));
            self.InvGamma_n = sparse(inv(Gamma_n));
            self.Ln = chol(self.InvGamma_n);%Store both invGamma and it's Cholesky
        end

        function [Hess, grad] = GetHessAndGrad(self, est, omega)
            if nargin < 3 || isempty(omega)
                omega = self.omega;
            end
            sigma = est.estimates{self.sigmaind} + 1i*omega*est.estimates{self.epsilonind};
            fres = self.SolveForwardVec(est, omega);
            if strcmp(self.mode, 'potential')
                diff = (fres - self.Iel);
            elseif strcmp(self.mode, 'current')
                diff = (fres - self.Uel);
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
            Hess = J'*self.InvGamma_n*J;%First order approximation of the forward problem
            grad = J'*self.InvGamma_n*diff;
            Hess = Estimate_Hess(Hess, est.GetSizes());
            grad = Estimate_vec(grad, est.GetSizes());

        end


        function sigma = PreProcessSigma(self, sigma)
            
            if length(self.sigmamin) == 1
                sigmin = [self.sigmamin; self.sigmamin];
            else
                sigmin = self.sigmamin;
            end
            if length(self.scales) == 2
                rsigma = self.scales(1)*real(sigma);
                isigma = self.scales(2)*imag(sigma);
                rsigma(rsigma<sigmin(1)) = sigmin(1);
                isigma(isigma<sigmin(2)) = sigmin(2);
                sigma = rsigma + 1i*isigma;
            elseif length(self.scales) > 1
                rsigma = self.scales(1:end/2).*real(sigma);
                isigma = self.scales(end/2+1:end).*imag(sigma);
                rsigma(rsigma<sigmin(1)) = sigmin(1);
                isigma(isigma<sigmin(2)) = sigmin(2);
                sigma = rsigma + 1i*isigma;
            else
                sigma = self.scales*sigma;
                sigma(real(sigma)<sigmin(1)) = sigmin(1) + 1i*imag(sigma(real(sigma)<sigmin(1)));
                sigma(imag(sigma)<sigmin(2)) = real(sigma(imag(sigma)<sigmin(2))) + 1i*sigmin(2);
            end
        
        end

    end

end