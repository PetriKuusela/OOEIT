classdef EITFEM < handle
%EITFEM class contains a FEM solver of the complete electrode model.
%Support of any kind of mesh is achieved by having the mesh in a separate
%class. More complex behaviour can be added by inheriting this class (for
%example complex EIT, MFEIT...)
%
%This class implements the functions required for it to be assigned to an
%inverse problem solver -class (e.g. SolverGN.m).
%
%Author: Petri Kuusela, 31.5.2024
    
    
    
    properties (SetObservable) %This block has to feature all the properties,
        %that affect the S and b parts of solving FEM. That way, we can
        %listen to changes and recalculate parts of the FEM matrix only when necessary.
        fmesh   %The mesh object for forward model
        C       %C matrix used to construct the basis function of electrode potentials/currents (nel x nel-1 array), nel = number of electrodes
        zeta    %Contact impedances of electrodes, nel x 1 array
        Iel     %The electrode currents (Iel and Uel can be either the measurements or the injection pattern, depending on 'mode')
        Uel     %The electrode potentials
        mode    %'current' or 'potential', which one is injected (the other ones are computed by solving CEM)
        Mpat    %Measurement pattern matrix. If measurements are potentials agains common ground, this is just 1
    end
    properties
        A       %The FEM matrix
        Ai      %The index decomposition used in caluclating the Jacobian
        Av
        intS    %A cell array containing integrals phi_i*phi_j along each electrode surfaces (one electrode per cell)
        intM    %An array (ng x nel) containing integral phi_i on surface of electrode j (ng = the number of nodes in the forward mesh)
        intB    %An array (ng x 1) containing the electrode areas in its elements
        mIncl   %A logical array (nel*ninj x 1) indicating which measurement values are used (ninj = number of injections in injection pattern)
        S       %The part of matrix A dependent on electrodes
        b       %Rhs of FEM equation A*Pot = b
        reCalc  %A flag for recalculating everything in the FEM matrix. This is automatically set to 1 if any of the (SetObservables)-block properties are modified.
        InvGamma%inverse of covariance matrix Gamma of the measurement values
        Ln      %chol(InvGamma)
        solVec  %The FEM solution A\b
        QC      %Matrix used to extract electrode values from Pot (i.e. Mpat*[0 C])
        sigmaMin%The minimum cutoff value for conductivity, all values below sigmamin will be set to sigmamin
        scales  %Scaling of each sigma value used before solving FEM. Useful if you want the estimate to have another scale than the actual conductivity does
        eps     %epsilon-correction epsilon, a vector added to every FEM result
        reCalcListener%A listener to listen if (SetObservables)-block parameters are modified
        zInd    %If estimate class (Estimate_vec) is used, which number is the contact impedances. If they are not estimated, zind = 0
        sigmaInd%If estimate class (Estimate_vec) is used, which number is the conductivity
        sigma   %Default conductivity, used if contact impedances are to be estimated without simultaneous conductivity estimation
        Iadded  %A flag to determine, if the user has remembered to add the measurement results before trying to solve the inverse problem
        Uadded  %A flag to determine, if the user has remembered to add the measurement results before trying to solve the inverse problem
    end %end properties
    
    
    methods
        
        function obj = EITFEM(fmesh)
            %Class constructor. Input (fmesh) is the forward mesh object 
            %containing all the mesh related computations. For example, it
            %can be a ForwardMesh1st object.
            obj.fmesh = fmesh;
            [obj.Ai, obj.Av] = fmesh.GradientMatrix();%These are used in the compution of the Jacobian
            [obj.intS, obj.intM, obj.intB] = fmesh.EITElectrodeTerms();%These integral values are used in computing the FEM matrix

            %populate necessary properties with default values
            %See the properties section for basic description of the
            %properties
            obj.C = [ones(1,fmesh.nEl-1);-eye(fmesh.nEl-1)];
            obj.sigmaMin = 1e-6;
            obj.Uel = eye(fmesh.nEl);
            obj.Uel = obj.Uel(:);%default injection pattern for potential injection
            obj.Iel = eye(fmesh.nEl) - [zeros(1,fmesh.nEl-1) 1; eye(fmesh.nEl-1) zeros(fmesh.nEl-1,1)];
            obj.Iel = obj.Iel(:);%default injection pattern for current injection
            obj.eps = 0;
            obj.mode = 'current';%current injection system as default
            obj.scales = 1;
            obj.zeta = 1e-6*ones(fmesh.nEl,1);
            obj.reCalc = 1;%On the first time, everything has to be calculated
            obj.Mpat = 1;
            obj.sigmaInd = 1;
            obj.zInd = 0;
            obj.sigma = ones(fmesh.ng, 1);

            %Set up a listener to listen for changes in the
            %(SetObservable)-block properties, to know when to set the
            %reCalc flag to 1.
            mc = metaclass(obj);
            metaprops = findobj([mc.Properties{:}], 'SetObservable', true);
            obj.reCalcListener = event.proplistener(obj, metaprops, 'PostSet', @obj.SetRecalc);
            obj.Iadded = 0;%Also keep track of whether measurement data has been added
            obj.Uadded = 0;
        end %end constructor
        
        function elval = SolveForward(self, sigma)
            %Solve the potentials or currents for given sigma. The output
            %is in matrix form, where each column corresponds to a single
            %injection.
                        
            %Multiply sigma by given scales and lift any values that are
            %too close to 0:
            sigma = self.PreProcessSigma(sigma);
                       
            %Start forming the EIT-matrix. This part is the part dependent
            %on the conductivity
            A0 = self.fmesh.SigmadPhiidPhij(sigma);
            
            %Get the injection pattern:
            if strcmp(self.mode, 'potential')
                inj = self.Uel;
            elseif strcmp(self.mode, 'current')
                inj = self.Iel;
            else
                error(['Unrecognized solver mode: ' self.mode]);
            end

            %Check if these values have already been calculated and have
            %not been modified since
            if self.reCalc
                %need to (re)calculate all values

                idr = repmat((1:self.fmesh.nEl)',1,self.fmesh.nEl-1);               %These indices are used just to put self.C into the correct 
                idc = repmat(self.fmesh.ng+(1:(self.fmesh.nEl-1)),self.fmesh.nEl,1);%position in the sparce matrix self.QC
                self.QC = self.Mpat'*sparse(idr,idc,self.C,self.fmesh.nEl,self.fmesh.ng+self.fmesh.nEl-1);

                %Compute integrals phi_i*phi_j/zeta along each electrode surfaces
                self.S = self.intS{1}/self.zeta(1);
                for ii=2:length(self.intS)
                    self.S = self.S + self.intS{ii}/self.zeta(ii);
                end

                %Compute the rest of the electrode dependent parts of the
                %FEM matrix
                if strcmp(self.mode, 'potential')
                    self.S = [self.S zeros(size(self.S,1), size(self.C,2)); -self.C'*diag(1./self.zeta)*self.intM' self.C'*self.C];
                    ninj = numel(inj)/self.fmesh.nEl;
                    self.b = zeros(size(self.S,1), ninj);
                    InjMat = reshape(inj, self.fmesh.nEl, ninj);
                    self.b(1:self.fmesh.ng,:) = self.intM*diag(1./self.zeta)*InjMat;
                    self.b(end-self.fmesh.nEl+2:end,:) = self.b(end-self.fmesh.nEl+2:end,:) - self.C'*(diag(self.intB./self.zeta)*InjMat);
                elseif strcmp(self.mode, 'current')
                    C2 = -self.C'*diag(1./self.zeta)*self.intM';
                    self.S = [self.S C2'; C2 +self.C'*diag(self.intB./self.zeta)*self.C];%check the sign of intB before using
                    ninj = numel(inj)/self.fmesh.nEl;
                    self.b = zeros(size(self.S,1), ninj);
                    InjMat = reshape(inj, self.fmesh.nEl, ninj);
                    self.b(end-self.fmesh.nEl+2:end,:) = self.C'*InjMat;
                else
                    error(['Unrecognized solver mode: ' self.mode]);
                end
                self.reCalc = 0;%This block does not have to be run again, until some relevant properties are changed
            end
            
            self.A = A0 + self.S; %Combine parts to make the full FEM matrix

            self.solVec = self.A\self.b;%Solve FEM

            elval = self.QC*self.solVec;%Extract the electrode values from the FEM solution (which contains also the potentials inside the domain)
            
        end %end solveForward

        function vec = SolveForwardVec(self, est)
            %Solve the FEM with given sigma.
            %output: vec = the currents or potentials in vector format
            %              computed using FEM
            %input:  est = the estimate. May be either a vector of
            %              conductivity values, or object of class
            %              EstimateVec, if multiple parameters are to be
            %              estimated.
            
            if isa(est, 'EstimateVec') %check if estimate-class is used
                if self.zInd > 0 && ~isempty(est.estimates{self.zInd})%Check if we want to estimate contact impedances
                    self.zeta = est.estimates{self.zInd};%If we estimate zeta, put the value from estimate-object on self.zeta
                end
                if isempty(est.estimates{self.sigmaInd})%Check if we want to estimate conductivity
                    vec = self.SolveForward(self.sigma);%Here, we do not want to estimate conductivity, so default value self.sigma is used
                else
                    vec = self.SolveForward(est.estimates{self.sigmaInd});%We get the conductivity from the estimate
                end
            else
                vec = self.SolveForward(est);%Estimate class is not used, so the estimate is just the conductivity
            end
            vec = vec(:);%Vectorize the output
            if ~isempty(self.mIncl)
                vec = vec(self.mIncl);%We can leave out some of the measurements by defining corresponding values at self.mIncl as false
            end
            vec = vec + self.eps;%Add a given constant (scalar or vector) to the results
            %This self.eps can be computed e.g. by SolveEpsilonCorrection.m
        end %end solveForwardVec
        
        function res = OptimizationFunction(self, sigma)
            %Calculate the residual of the forward problem, which is to be
            %minimized (in addition to the regularization) when solving the
            %inverse problem. This function is called by the inverse
            %problem solver (e.g. SolverGN.m)        
           
            elVal = self.SolveForwardVec(sigma);%The electrode potentials or currents
            if strcmp(self.mode, 'potential')
                if ~self.Iadded%Check if measurements have not been added, i.e. we try to solve with the default injection pattern as measurements
                    warning('Attempting optimization with the default injection pattern as data! Did you forget to load the measurement data to EITFEM object?');
                    self.Iadded = 1;%suppress warning after first time
                end

                %This is the residual value to be minimized
                res = 0.5*sum((self.Ln*(abs(self.Iel - elVal))).^2);

            elseif strcmp(self.mode, 'current')
                if ~self.Uadded%Check if measurements have not been added, i.e. we try to solve with the default injection pattern as measurements
                    warning('Attempting optimization with the default injection pattern as data! Did you forget to load the measurement data to EITFEM object?');
                    self.Uadded = 1;%suppress warning after first time
                end

                %This is the residual value to be minimized
                res = 0.5*sum((self.Ln*(abs(self.Uel - elVal))).^2);
                
            else
                error(['Unrecognized solver mode: ' self.mode]);
            end
        end
        
        function SetInvGamma(self, constNoise, relNoise)
            %Calculate and set the inverse of covariance matrix based on
            %the noise levels given as arguments.
            %Input: const_noise = a constant noise level for all
            %       measurements. This noise is scaled to the
            %       difference of minimum and maximum measurement values
            %
            %       rel_noise(optional) = relative noise level of the
            %       measurements. This is multiplied by the absolute value
            %       of each measurement to get the noise level of that
            %       measurement. 
            %
            %       The two types of noise above are added together.
            %NOTE: This does not add noise to the solution of the forward
            %problem! This just calculates the Weighing matrix used in the
            %inverse problem.
            
            %Check the optional arguments:
            if nargin < 3 || isempty(relNoise)
                relNoise = 0;
            end
            
            %Determine which (current or potential) is the dependent
            %variable
            if strcmp(self.mode, 'potential')
                dMeas = self.Iel;%dependent measurements
            elseif strcmp(self.mode, 'current')
                dMeas = self.Uel;%dependent measurements
            else
                error(['Unrecognized solver mode: ' self.mode]);
            end    
            
            %Calculating the model variance of the noise
            %Constant noise for all measurements:
            varMeas = (constNoise*(max(dMeas)-min(dMeas)))^2;
            %Add noise relative to the abs of measurement value:
            varMeas = varMeas + (relNoise*(abs(dMeas))).^2;
            %Assume no cross-correlations, hence a diagonal covariance mat:
            Gamma = diag(varMeas(:));
            self.InvGamma = sparse(inv(Gamma));%This is the data precision matrix
            self.Ln = chol(self.InvGamma);%Store both invGamma and it's Cholesky
        
        end % end SetInvGamma
        
        function [Hess, grad] = GetHessAndGrad(self, est)
            %Calculate the Hess-matrix and gradient used in solving the
            %inverse problem. This function is called by the inverse
            %problem solver (e.g. SolverGN.m)
            %
            %The Hess-matrix returned by this function is not the full
            %Hess-matrix, but the one used in Gauss-Newton approximation,
            %i.e. only 1st derivatives of the CEM are taken into account.
            %
            %input: est =  the estimate. May be either a vector of
            %              conductivity values, or object of class
            %              EstimateVec, if multiple parameters are to be
            %              estimated.
            
            %solve the electrode values
            elVal = self.SolveForwardVec(est);
            %and find the difference to the measurements
            if strcmp(self.mode, 'potential')
                diff = (elVal - self.Iel);
            elseif strcmp(self.mode, 'current')
                diff = (elVal - self.Uel);
            else
                error(['Unrecognized solver mode: ' self.mode]);
            end


            if isa(est, 'EstimateVec') %Check if estimate object is used
                HC = cell(length(est.estimates));%Initialize the Hess and gradient cell arrays for the return objects
                gC = cell(length(est.estimates),1);
                if ~isempty(est.estimates{self.sigmaInd})%Check are we estimating conductivity
                    Js = self.Jacobian(est, 1);%compute the Jacobian
                    Js = Js*diag(self.scales);%The scaling done before solving the FEM (in self.PreProcessSigma()) has to be taken into account here.
                    HC{self.sigmaInd, self.sigmaInd} = Js'*self.InvGamma*Js;%First order approximation of the forward problem
                    gC{self.sigmaInd} = Js'*self.InvGamma*diff;%gradient
                end
                if self.zInd > 0 && ~isempty(est.estimates{self.zInd})%Check are we estimating the contact impedances
                    Jz = self.JacobianZ(est, 1);%The Jacobian of the contact impedances
                    HC{self.zInd, self.zInd} = Jz'*self.InvGamma*Jz;%First order approximation of the forward problem
                    gC{self.zInd} = Jz'*self.InvGamma*diff;%gradient
                end
                if self.zInd > 0 && ~isempty(est.estimates{self.sigmaInd}) && ~isempty(est.estimates{self.zInd})
                    %estimating both contact impedance and conductivity, so
                    %need the cross terms as well
                    HC{self.sigmaInd, self.zInd} = Js'*self.InvGamma*Jz;
                end
                Hess = EstimateHess(HC);%Create the output objects
                grad = EstimateVec(gC);
            else %Estimate objects are not used, the estimate is just a conductivity vector
                J = self.Jacobian(est, 1);
                J = J*diag(self.scales);%The scaling done before solving the FEM (in self.PreProcessSigma()) has to be taken into account here.
                Hess = J'*self.InvGamma*J;%First order approximation of the forward problem
                grad = J'*self.InvGamma*diff;
            end
            
        end
        
        function J = Jacobian(self, sigma, alreadyComputed)
            %Calculate the Jacobian of the forward problem at est (i.e.
            %sigma). The Jacobian = the derivatives of the measurements
            %(electrode currents or potentials) w.r.t. the internal
            %conductivity values.

            %alreadyComputed is used to determine whether we need to
            %compute the forward problem again, to make sure self.solVec
            %corresponds to the conductivity given. Defaul is to recompute
            %to make sure.
            if nargin < 3 || isempty(alreadyComputed)
                alreadyComputed = 0;
            end

            if ~alreadyComputed%Have to compute the forward problem to make sure self.solVec is up-to-date
                self.SolveForward(sigma);
            end

            ng = length(self.Ai);%number of nodes in the forward mesh
            c = size(self.QC,1); %number of electrodes
            d = size(self.solVec,2); %number of injections

            %The Jacobian is essentially computed as
            %-self.QC*inv(self.A)*dA*self.solVec, where dA is the derivative
            %of the FEM matrix A w.r.t. sigma. Since dA has only a few 
            %indices with non-zero entries, for computational
            %performance, the calculations have been split in the following
            %way.

            Jleft = -self.QC/self.A;%These parts have to be computed only once
            Jright = self.solVec;

            %Initialize the jacobian matrix
            Js = zeros(c*d,ng);

            for ii=1:ng
              Jid = self.Ai{ii};%The non-zero indices of dA/dsigma for node ii
              Jtemp   = Jleft(:,Jid)*self.Av{ii}*Jright(Jid,:);%This corresponds to a single column of the Jacobian
              Js(:,ii) = Jtemp(:);
            end
            
            if ~isempty(self.mIncl)%Check if we want to leave some measurements out
                Js = Js(self.mIncl,:);
            end
            J = self.fmesh.JacobianFtoI(Js);%If a separate mesh is used for conductivity, transform the Jacobian from the FEM basis to condcutivity basis
            
        end

        function Js = JacobianZ(self, est, alreadyComputed)
            % Computes the contact impedance Jacobian J = d(measurements) / d(zeta)

            %alreadyComputed is used to determine whether we need to
            %compute the forward problem again, to make sure self.Pot
            %corresponds to the conductivity given. Defaul is to recompute
            %to make sure.
            if nargin < 3 || isempty(alreadyComputed)
                alreadyComputed = 0;
            end

            if ~alreadyComputed
                self.SolveForwardVec(est);
            end

            ng = length(self.zeta);%number of nodes in the forward mesh
            c = size(self.QC,1);   %number of electrodes
            d = size(self.solVec,2);  %number of injections

            if strcmp(self.mode, 'potential')%In the case of potential injection, we need the injection matrix
                Inj = reshape(self.Uel, c, d);
            end

            %For more comments on computing the Jacobian, see the basic
            %Jacobian above
            Jleft  = self.QC/self.A;
            Jright = self.solVec;

            Js = zeros(c*d,ng);

            for ii=1:ng
                rzeta = zeros(size(self.zeta));
                rzeta(ii) = 1;
                tempM = self.intM*diag(rzeta);
                if strcmp(self.mode, 'potential')
                    tempS = [self.intS{ii} zeros(size(self.intS{ii},1), size(self.C,2));...
                                    -self.C'*tempM' zeros(size(self.C,2))];
                    ninj = numel(Inj)/(self.fmesh.nEl);
                    tb = zeros(size(tempS,1), ninj);
                    tb(1:self.fmesh.ng,:) = tempM*Inj;
                    tb(end-self.fmesh.nEl+2:end,:) = ...
                                              - self.C'*(diag(self.intB.*rzeta)*Inj);
                    Jtemp   = -1/est.estimates{self.zInd}(ii)^2*Jleft*(tb - tempS*Jright);
                elseif strcmp(self.mode, 'current')
                    tempB = zeros(size(self.intB));
                    tempB(ii) = self.intB(ii);
                    tempS = [self.intS{ii} -tempM*self.C;...
                                -self.C'*tempM' self.C'*diag(tempB)*self.C];
                    Jtemp   = 1/est.estimates{self.zInd}(ii)^2*Jleft*tempS*Jright;
                else
                    error(['Unrecognized solver mode: ' self.mode]);
                end
                Js(:,ii) = Jtemp(:);
            end

            if ~isempty(self.mIncl)
                Js = Js(self.mIncl,:);
            end

        end
        
        function sigma = PreProcessSigma(self, sigma)
            
            sigma = sigma.*self.scales;%Apply given scales to sigma
            sigma(sigma<self.sigmaMin) = self.sigmaMin;%Lift all the negative (or almost) values to the minimum value
        
        end
        
        function Plot(self, sigma)
            %Used to plot how the FEM results fit with the measurement data
            
            elVal = self.SolveForwardVec(sigma);
            if strcmp(self.mode, 'potential')
                plot(1:length(elVal), elVal, 'r-', 1:length(self.Iel), self.Iel, 'b-');
            elseif strcmp(self.mode, 'current')
                plot(1:length(elVal), elVal, 'r-', 1:length(self.Uel), self.Uel, 'b-');
            else
                error(['Unrecognized solver mode: ' self.mode]);
            end
            legend('Forward', 'Measurement');
        end
       
        function SetRecalc(obj, ~, ~)
            %This function is used by a listener to set up recalc-flag
            %whenever all the elements of the FEM matrix need to be
            %recalculated.
            obj.reCalc = 1;
        end

        function set.Iel(self, Iel)
            %Make note when Iel is set, so that we know it is not anymore
            %just the default injection pattern.
            %These could possibly cause trouble whever loading the solver from a file
            self.Iel = Iel;
            self.Iadded = 1;
            if isempty(self.InvGamma) && strcmp(self.mode, 'potential')
                self.SetInvGamma(1e-3, 3e-2)%Also compute some form of data-precision matrix
            end
        end

        function set.Uel(self, Uel)
            %Make note when Uel is set, so that we know it is not anymore
            %just the default injection pattern
            %These could possibly cause trouble whever loading the solver from a file
            self.Uel = Uel;
            self.Uadded = 1;
            if isempty(self.InvGamma) && strcmp(self.mode, 'current')
                self.SetInvGamma(1e-3, 3e-2)%Also compute some form of data-precision matrix
            end
        end

    end

    
end