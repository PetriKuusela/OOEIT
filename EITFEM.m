classdef EITFEM < handle
%EITFem class contains a FEM solver of the complete electrode model.
%Support of any kind of mesh is achieved by having the mesh in a separate
%class. More complex behaviour can be added by inheriting this class (for
%example complex EIT, MFEIT...)
%
%This class implements the functions required for it to be assigned to an
%inverse problem solver -class (e.g. SolverLinesearch.m).
%
%Author: Petri Kuusela, 5.4.2024
    
    
    
    properties (SetObservable) %This block has to feature all the properties, that affect the S and b parts of solving FEM.
        fmesh   %forward mesh
        C       %C matrix used to construct the basis function of electrode potentials/currents (nel x nel-1 array)
        zeta    %Contact impedances of electrodes, nel x 1 array
        Iel     %The electrode currents
        Uel     %The electrode potentials
        mode    %'current' or 'potential', which one is injected (the other one is solved by solving CEM)
        Mpat    %Measurement pattern matrix. If measurements are potentials agains common ground, this is just 1
    end
    properties
        A       %FEM matrix
        Ai      %The index decomposition used in caluclating the Jacobian
        Av
        intS    %integral values used for the FEM matrix
        intM    %integral values used for the FEM matrix
        intB    %integral values used for the FEM matrix
        vincl   %A logical array (ng x 1) indicating which measurement values are used
        S       %The part of matrix A dependent on electrodes
        b       %Rhs of FEM equation A*Pot = b
        recalc  %A flag for recalculating everything in the FEM matrix. e.g. set to 1 if zetas are modified.
        InvGamma_n %inverse of covariance matrix Gamma of the measurement values
        Ln      %chol(InvGamma_n)
        Pot     %The FEM solution A\b
        QC      %Matrix used to extract solutions from Pot (i.e. Mpat*[0 0; 0 C])
        sigmamin%The minimum cutoff value for conductivity, all values below sigmamin will be set to sigmamin
        scales  %Scaling of each sigma value used before solving FEM. If len(scales) == 2 the first value scales real(sigma) and second imag(sigma)
        eps     %epsilon-correction epsilon, added to every FEM result
        recalclistener
        zind
        sigmaind
        sigma
    end %end properties
    
    
    methods
        
        function obj = EITFEM(fmesh)
            %Class constructor. Input (fmesh) is the forward mesh of
            %the FEM solver (see properties-section for details).
            obj.fmesh = fmesh;
            [obj.Ai, obj.Av] = fmesh.GradientMatrix();
            
            %populate necessary properties with default values
            [obj.intS, obj.intM, obj.intB] = fmesh.EITElectrodeTerms();
            obj.C = [ones(1,fmesh.nel-1);-eye(fmesh.nel-1)];
            obj.sigmamin = 1e-6;
            obj.vincl = true(fmesh.nel*fmesh.nel,1);
            obj.mode = 'potential';%potential injection system
            obj.Uel = diag(ones(fmesh.nel,1));
            obj.Uel = obj.Uel(:);
            obj.eps = 0;
            obj.scales = 1;
            obj.zeta = 1e-9*ones(fmesh.nel,1);
            obj.recalc = 1;
            obj.Mpat = 1;
            mc = metaclass(obj);
            metaprops = findobj([mc.Properties{:}], 'SetObservable', true);
            obj.recalclistener = event.proplistener(obj, metaprops, 'PostSet', @obj.SetRecalc);
            obj.sigmaind = 1;
            obj.zind = 0;
            obj.sigma = ones(fmesh.ng, 1);
        end %end constructor
        
        function elval = SolveForward(self, sigma)
            %Solve the potentials or currents for given sigma
                        
            %Multiply sigma by given scales and lift any values that are
            %too close to 0:
            sigma = self.PreProcessSigma(sigma);
                       
            %Start forming the EIT-matrix
            A0 = self.fmesh.SigmadPhiidPhij(sigma);
            
            if strcmp(self.mode, 'potential')
                Inj = self.Uel;
            elseif strcmp(self.mode, 'current')
                Inj = self.Iel;
            else
                error(['Unrecognized solver mode: ' self.mode]);
            end

            if self.recalc
                idr = repmat((1:self.fmesh.nel)',1,self.fmesh.nel-1);
                idc = repmat(self.fmesh.ng+(1:(self.fmesh.nel-1)),self.fmesh.nel,1);
                self.QC = self.Mpat'*sparse(idr,idc,self.C,self.fmesh.nel,self.fmesh.ng+self.fmesh.nel-1);

                self.S = self.intS{1}/self.zeta(1);
                for ii=2:length(self.intS)
                    self.S = self.S + self.intS{ii}/self.zeta(ii);
                end
                if strcmp(self.mode, 'potential')
                    self.S = [self.S zeros(size(self.S,1), size(self.C,2)); -self.C'*diag(1./self.zeta)*self.intM' self.C'*self.C];
                    ninj = numel(Inj)/self.fmesh.nel;
                    self.b = zeros(size(self.S,1), ninj);
                    InjMat = reshape(Inj, self.fmesh.nel, ninj);
                    self.b(1:self.fmesh.ng,:) = self.intM*diag(1./self.zeta)*InjMat;
                    self.b(end-self.fmesh.nel+2:end,:) = self.b(end-self.fmesh.nel+2:end,:) - self.C'*(diag(self.intB./self.zeta)*InjMat);
                elseif strcmp(self.mode, 'current')
                    C2 = -self.C'*diag(1./self.zeta)*self.intM';
                    self.S = [self.S C2'; C2 +self.C'*diag(self.intB./self.zeta)*self.C];%check the sign of intB before using
                    ninj = numel(Inj)/self.fmesh.nel;
                    self.b = zeros(size(self.S,1), ninj);
                    InjMat = reshape(Inj, self.fmesh.nel, ninj);
                    self.b(end-self.fmesh.nel+2:end,:) = self.C'*InjMat;
                else
                    error(['Unrecognized solver mode: ' self.mode]);
                end
                self.recalc = 0;
            end
            
            self.A = A0 + self.S;

            self.Pot = self.A\self.b;

            elval = self.QC*self.Pot;
            
        end %end solveForward

        function vec = SolveForwardVec(self, est)
            %Solve FEM with given sigma.
            %output: vec = the currents or potentials in vector format
            %              computed using FEM
            
            if isa(est, 'Estimate_vec')
                if self.zind > 0 && ~isempty(est.estimates{self.zind})
                    self.zeta = est.estimates{self.zind};
                end
                if isempty(est.estimates{self.sigmaind})
                    vec = self.SolveForward(self.sigma);
                else
                    vec = self.SolveForward(est.estimates{self.sigmaind});
                end
            else
                vec = self.SolveForward(est);
            end
            vec = vec(:);
            vec = vec(self.vincl);
            vec = vec + self.eps;%Add a given constant (scalar or vector) to the results
        end %end solveForwardVec
        
        function res = OptimizationFunction(self, sigma)
            %Calculate the residual of the forward problem, which is to be
            %minimized (in addition to the regularization) when solving the
            %inverse problem. This function is called by the inverse
            %problem solver (e.g. SolverLinesearch.m)        
           
            elval = self.SolveForwardVec(sigma);
            if strcmp(self.mode, 'potential')
                res = 0.5*sum((self.Ln*(abs(self.Iel - elval))).^2);
            elseif strcmp(self.mode, 'current')
                res = 0.5*sum((self.Ln*(abs(self.Uel - elval))).^2);
            else
                error(['Unrecognized solver mode: ' self.mode]);
            end
        end
        
        function SetInvGamma(self, meas_noise_coef_e, meas_noise_coef2)
            %Calculate and set the inverse of covariance matrix based on
            %the noise levels given as arguments.
            %Input: meas_noise_coef_e = a constant noise level for all
            %       measurements. This coefficient is scaled to the
            %       difference of minimum and maximum measurement values to
            %       get the actual noise level.
            %
            %       meas_noise_coef2(optional) = relative noise level of the
            %       measurements. This is multiplied by the absolute value
            %       of each measurement to get the noise level of that
            %       measurement. 
            %
            %       The two types of noise above are added together.
            %NOTE: This does not add noise to the solution of the forward
            %problem! This just calculates the Weighing matrix used in the
            %inverse problem.
            
            %Check the optional arguments:
            if nargin < 3 || isempty(meas_noise_coef2)
                meas_noise_coef2 = 0;
            end
            
            if strcmp(self.mode, 'potential')
                Dmeas = self.Iel;
            elseif strcmp(self.mode, 'current')
                Dmeas = self.Uel;
            else
                error(['Unrecognized solver mode: ' self.mode]);
            end    
            
            %Calculating the model variance of the noise
            %Constant noise for all measurements:
            var_meas = (meas_noise_coef_e*(max(max(Dmeas))-min(min(Dmeas))))^2;
            %Noise level relative to the abs of measurement value:
            var_meas = var_meas + (meas_noise_coef2*(abs(Dmeas) - min(min(abs(Dmeas))))).^2;
            %Assume no cross-correlations, hence a diagonal covariance mat:
            Gamma_n = diag(var_meas(:));
            self.InvGamma_n = sparse(inv(Gamma_n));
            self.Ln = chol(self.InvGamma_n);%Store both invGamma and it's Cholesky
        
        end % end SetInvGamma
        
        function [Hess, grad] = GetHessAndGrad(self, est)
            %Calculate the Hess-matrix and gradient used in solving the
            %inverse problem. This function is called by the inverse
            %problem solver (e.g. SolverLinesearch.m)
            %
            %The Hess-matrix returned by this function is not the full
            %Hess-matrix, but the one used in Gauss-Newton approximation,
            %i.e. only 1st derivatives of the CEM are taken into account.
            
            fres = self.SolveForwardVec(est);
            if strcmp(self.mode, 'potential')
                diff = (fres - self.Iel);
            elseif strcmp(self.mode, 'current')
                diff = (fres - self.Uel);
            else
                error(['Unrecognized solver mode: ' self.mode]);
            end
            if isa(est, 'Estimate_vec')
                HC = cell(length(est.estimates));
                gC = cell(length(est.estimates),1);
                if ~isempty(est.estimates{self.sigmaind})
                    Js = self.Jacobian(est, 1);
                    Js = Js*diag(self.scales);%The scaling done before solving the FEM has to be taken into account here.
                    HC{self.sigmaind, self.sigmaind} = Js'*self.InvGamma_n*Js;%First order approximation of the forward problem
                    gC{self.sigmaind} = Js'*self.InvGamma_n*diff;
                end
                if self.zind > 0 && ~isempty(est.estimates{self.zind})
                    Jz = self.Jacobian_z(est, 1);
                    HC{self.zind, self.zind} = Jz'*self.InvGamma_n*Jz;%First order approximation of the forward problem
                    gC{self.zind} = Jz'*self.InvGamma_n*diff;
                end
                if self.zind > 0 && ~isempty(est.estimates{self.sigmaind}) && ~isempty(est.estimates{self.zind})
                    HC{self.sigmaind, self.zind} = Js'*self.InvGamma_n*Jz;
                end
                Hess = Estimate_Hess(HC);
                grad = Estimate_vec(gC);
            else
                J = self.Jacobian(est, 1);
                J = J*diag(self.scales);%The scaling done before solving the FEM has to be taken into account here.
                Hess = J'*self.InvGamma_n*J;%First order approximation of the forward problem
                grad = J'*self.InvGamma_n*diff;
            end
            
        end
        
        function J = Jacobian(self, est, alreadyComputed)
            %Calculate the Jacobian of the forward problem at est (i.e.
            %sigma)
            if nargin < 3 || isempty(alreadyComputed)
                alreadyComputed = 0;
            end

            b = length(self.Ai);
            c = size(self.QC,1); 
            d = size(self.Pot,2); 

            if ~alreadyComputed
                self.SolveForward(est);
            end

            Jleft = -self.QC/self.A;

            Jright = self.Pot;

            Js = zeros(c*d,b);

            for ii=1:b
              Jid = self.Ai{ii};

              Jtemp   = Jleft(:,Jid)*self.Av{ii}*Jright(Jid,:);
              Js(:,ii) = Jtemp(:);
            end
            
            Js = Js(self.vincl,:);
            J = self.fmesh.JacobianFtoI(Js);
            
        end

        function Js = Jacobian_z(self, est, alreadyComputed)
            % Computes the Jacobian J = d(measurements) / d(zeta)
            if nargin < 3 || isempty(alreadyComputed)
                alreadyComputed = 0;
            end

            b = length(self.zeta);
            c = size(self.QC,1); 
            d = size(self.Pot,2); 

            if ~alreadyComputed
                self.SolveForward(est);
            end

            if strcmp(self.mode, 'potential')
                Inj = reshape(self.Uel, c, d);
            end
            Jleft  = self.QC/self.A;

            Jright = self.Pot;

            Js = zeros(c*d,b);

            for ii=1:b
                rzeta = zeros(size(self.zeta));
                rzeta(ii) = 1;
                tempM = self.intM*diag(rzeta);
                if strcmp(self.mode, 'potential')
                    tempS = [self.intS{ii} zeros(size(self.intS{ii},1), size(self.C,2));...
                                    -self.C'*tempM' zeros(size(self.C,2))];
                    ninj = numel(Inj)/(self.fmesh.nel);
                    tb = zeros(size(tempS,1), ninj);
                    tb(1:self.fmesh.ng,:) = tempM*Inj;
                    tb(end-self.fmesh.nel+2:end,:) = ...
                                              - self.C'*(diag(self.intB.*rzeta)*Inj);
                    Jtemp   = -1/est.estimates{self.zind}(ii)^2*Jleft*(tb - tempS*Jright);
                elseif strcmp(self.mode, 'current')
                    tempB = zeros(size(self.intB));
                    tempB(ii) = self.intB(ii);
                    tempS = [self.intS{ii} -tempM*self.C;...
                                -self.C'*tempM' self.C'*diag(tempB)*self.C];
                    Jtemp   = 1/est.estimates{self.zind}(ii)^2*Jleft*tempS*Jright;
                else
                    error(['Unrecognized solver mode: ' self.mode]);
                end
                Js(:,ii) = Jtemp(:);
            end

        end
        
        function sigma = PreProcessSigma(self, sigma)
            
            sigma = sigma.*self.scales;
            sigma(sigma<self.sigmamin) = self.sigmamin;
        
        end
        
        function Plot(self, sigma)
            %Used to plot how the FEM results fit with the measurement data
            
            elval = self.SolveForwardVec(sigma);
            if strcmp(self.mode, 'potential')
                plot(1:length(elval), elval, 'r-', 1:length(self.Iel), self.Iel, 'b-');
            elseif strcmp(self.mode, 'current')
                plot(1:length(elval), elval, 'r-', 1:length(self.Uel), self.Uel, 'b-');
            else
                error(['Unrecognized solver mode: ' self.mode]);
            end
            legend('Forward', 'Measurement');
        end
       
        function SetRecalc(obj, ~, ~)
            %This function is used by a listener to set up recalc-flag
            %whenever all the elements of the FEM matrix need to be
            %recalculated.
            obj.recalc = 1;
        end

    end

    
end