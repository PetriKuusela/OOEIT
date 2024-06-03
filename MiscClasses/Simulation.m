function [Umeas, Imeas, Umeas_i, Imeas_i, e] = Simulation(fm, sigma_i, sigma, z, fsolver, mode, err)
%A simple simulation that will return synthetic data. meas_i refer to
%initial reference measurements, whereas Umeas and Imeas are the situation
%we are actually interested in. if sigma_i = [], reference measurements are
%not calculated.
%
%Author: Petri Kuusela 5.4.2024

%Use default values for arguments that are missing:
if nargin < 3 || isempty(sigma)
    sigma = ones(fm.nginv,1);
end
if nargin < 4 || isempty(z)
    z = 1e-6*ones(fm.nel,1);
end
if nargin < 6 || isempty(mode)
    if nargin > 4 && ~isempty(fsolver)
        mode = fsolver.mode;
    else
        mode = 'current';
    end
end
if nargin < 7 || isempty(err)
    err = [1e-6 1e-3 1e-4 1e-2];
end

    %These are the noise parameters
    noiserel = err(1);
    noiseabs = err(2);
    esystematicrel = err(3);
    esystematicabs = err(4);

    if isempty(fsolver)
        fsolver = EITFEM(fm);
        fsolver.zeta = z;
        fsolver.mode = mode;
        fsolver.sigmamin = 1e-9;
        if strcmp(mode, 'potential')
            fsolver.Uel = eye(length(fm.E));%Injection pattern
            fsolver.Uel= fsolver.Uel(:);
        elseif strcmp(mode, 'current')
            Imeas = eye(length(fm.E));%Create the injection pattern
            Imeas(2:end,1:end-1) = Imeas(2:end,1:end-1) - eye(length(fm.E)-1);
            Imeas(1,end) = -1;
            fsolver.Iel = Imeas(:);
            fsolver.Iel = fsolver.Iel*1e-3;%Injected currents are usually order of mA
        end    
    end

    if ~isempty(sigma_i)
        %reference measurements:
        if length(sigma_i) < 2
            sigma_i = [sigma_i;sigma_i];
        end

        if strcmp(mode, 'potential')
            Imeas_i = fsolver.SolveForwardVec(sigma_i);%these are the results
            Umeas_i = fsolver.Uel;
            %add noise and error:
            esys = randn(length(Imeas_i),1)*esystematicabs*(max(Imeas_i)-min(Imeas_i));
            esysrel = randn(length(Imeas_i),1)*esystematicrel;
            Imeas_i = Imeas_i.*(1+esysrel) + esys;
            Imeas_i = Imeas_i.*(1+randn(length(Imeas_i),1)*noiserel) + randn(length(Imeas_i),1)*noiseabs*(max(Imeas_i)-min(Imeas_i));
        elseif strcmp(mode, 'current')
            Umeas_i = fsolver.SolveForwardVec(sigma_i);%these are the results
            Imeas_i = fsolver.Iel;
            %add noise and error:
            esys = randn(length(Umeas_i),1)*esystematicabs*(max(Umeas_i)-min(Umeas_i));
            esysrel = randn(length(Umeas_i),1)*esystematicrel;
            Umeas_i = Umeas_i.*(1+esysrel) + esys;
            Umeas_i = Umeas_i.*(1+randn(length(Umeas_i),1)*noiserel) + randn(length(Umeas_i),1)*noiseabs*(max(Umeas_i)-min(Umeas_i));
        end

    end %end reference measurements
        
        
    %The actual measurements:
    if strcmp(mode, 'potential')
        Imeas = fsolver.SolveForwardVec(sigma);%these are the results
        Umeas = fsolver.Uel;
        %add noise and error:
        if ~exist('esys')%there have been no homogeneous measurements where these have already been calculated
            esys = randn(length(Imeas),1)*esystematicabs*(max(Imeas)-min(Imeas));
            esysrel = randn(length(Imeas),1)*esystematicrel;
        end%If previous esys and esysrel exist use them
        Imeas = Imeas.*(1+esysrel) + esys;
        Imeas = Imeas.*(1+randn(length(Imeas),1)*noiserel) + randn(length(Imeas),1)*noiseabs*(max(Imeas)-min(Imeas));
    elseif strcmp(mode, 'current')
        Umeas = fsolver.SolveForwardVec(sigma);%these are the results
        Imeas = fsolver.Iel;
        %add noise and error:
        if ~exist('esys')%there have been no homogeneous measurements where these have already been calculated
            esys = randn(length(Umeas),1)*esystematicabs*(max(Umeas)-min(Umeas));
            esysrel = randn(length(Umeas),1)*esystematicrel;
        end%If previous esys and esysrel exist use them
        Umeas = Umeas.*(1+esysrel) + esys;
        Umeas = Umeas.*(1+randn(length(Umeas),1)*noiserel) + randn(length(Umeas),1)*noiseabs*(max(Umeas)-min(Umeas));
    end
    
    %return also the systematic error:
    e = [esys esysrel];
    if ~exist('Umeas_i')%reference measurements were not done
        Umeas_i = 0;
        Imeas_i = 0;
    end
    
    
end