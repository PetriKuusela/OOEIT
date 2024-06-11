classdef PriorPositivityParabolic < handle
    %PositivityParabolic class contains a simple parabolic positivity
    %constraint and the required functions to be called by an inverse
    %problems solver (e.g. SolverLinesearch.m).
    %
    %Author: Petri Kuusela 5.4.2024
    properties
        a %the constraint function value will be a*(x-minAcceptable)^2 for all x < minAcceptable
        minAcceptable%The minimum value of estimate that is not penalized
        omitInd
    end
    methods
        function obj = PriorPositivityParabolic(minAcceptable, valAt0)
            %Class constructor.
            %Input: minAcceptable = The minimum value of estimate which is
            %                       not penalized.
            %       valAt0 = What is the functional value at 0
            obj.minAcceptable = minAcceptable;
            obj.a = valAt0./(minAcceptable.^2);
        end
        function res = OptimizationFunction(self, est)
            if isa(est, 'EstimateVec')
                res = 0;
                for ii = 1:length(est.estimates)
                    if ismember(ii, self.omitInd)
                        continue;
                    end
                    if length(self.a) < ii
                        self.a(ii) = self.a(1);
                    end
                    if length(self.minAcceptable) < ii
                        self.minAcceptable(ii) = self.minAcceptable(1);
                    end
                    res = res + self.a(ii)*sum((est.estimates{ii}(est.estimates{ii}<self.minAcceptable(ii)) - self.minAcceptable(ii)).^2);
                end
            else
                res = self.a(1)*sum((est(est<self.minAcceptable(1)) - self.minAcceptable(1)).^2);
            end
        end
        function [Hess, grad] = GetHessAndGrad(self, est)

            if isa(est, 'EstimateVec')
                HC = cell(length(est.estimates));
                gC = cell(length(est.estimates),1);
                for ii = 1:length(est.estimates)
                    if ismember(ii, self.omitInd)
                        continue;
                    end
                    if length(self.a) < ii
                        self.a(ii) = self.a(1);
                    end
                    if length(self.minAcceptable) < ii
                        self.minAcceptable(ii) = self.minAcceptable(1);
                    end
                    gC{ii} = zeros(length(est.estimates{ii}),1);
                    select = est.estimates{ii}<self.minAcceptable(ii);
                    gC{ii}(select) = 2*self.a(ii)*(est.estimates{ii}(select) - self.minAcceptable(ii));
                    Hessvec = zeros(length(est.estimates{ii}),1);
                    Hessvec(select) = 2*self.a(ii);
                    HC{ii,ii} = diag(Hessvec);
                end
                Hess = EstimateHess(HC);
                grad = EstimateVec(gC);
            else
                grad = zeros(length(est),1);
                select = est<self.minAcceptable;
                grad(select) = 2*self.a*(est(select) - self.minAcceptable);
                Hessvec = zeros(length(est),1);
                Hessvec(select) = 2*self.a;
                Hess = diag(Hessvec);
            end

        end
        
    end
end
        