classdef EstimateHess
%This is a class used for estimating multiple parameters simultaneously,
%e.g. conductivity + contact impedances. Used in conjunction with
%EstimateVec.
%
%This class has implementations for the basic operations (like + and \)
%required for it to function in an inverse solver, e.g. SolverGN, just like
%a matrix would.

    properties
        estimates
    end


    methods

        function obj = EstimateHess(estimates, sizes)
            if nargin == 1
                obj.estimates = estimates;
            elseif nargin == 2
                est = cell(length(sizes));
                libound = 1;
                ljbound = 1;
                for ii = 1:length(sizes)
                    for ij = ii:length(sizes)
                        est{ii,ij} = estimates(libound:libound+sizes(ii)-1, ljbound:ljbound+sizes(ij)-1);
                        ljbound = ljbound + sizes(ij);
                    end
                    libound = libound + sizes(ii);
                    ljbound = sum(sizes(1:ii))+1;
                end
                obj.estimates = est;
            end
        end

        function r = plus(o1, o2)
            if isa(o1, 'EstimateHess') && isa(o2, 'EstimateHess')
                newest = cell(size(o1.estimates,1));
                for ii = 1:size(o1.estimates,1)
                    for ij = 1:size(o1.estimates,2)
                        if isempty(o1.estimates{ii,ij})
                            newest{ii,ij} = o2.estimates{ii,ij};
                        elseif isempty(o2.estimates{ii,ij})
                            newest{ii,ij} = o1.estimates{ii,ij};
                        else
                            newest{ii,ij} = o1.estimates{ii,ij} + o2.estimates{ii,ij};
                        end
                    end
                end
            elseif isa(o1, 'EstimateHess')
                newest = cell(size(o1.estimates,1));
                for ii = 1:size(o1.estimates,1)
                    for ij = 1:size(o1.estimates,2)
                        newest{ii,ij} = o1.estimates{ii,ij} + o2;
                    end
                end
            else
                newest = cell(size(o2.estimates,1));
                for ii = 1:size(o2.estimates,1)
                    for ij = 1:size(o2.estimates,2)
                        newest{ii,ij} = o2.estimates{ii,ij} + o1;
                    end
                end
            end
            r = EstimateHess(newest);
        end

        function r = minus(o1, o2)
            if isa(o1, 'EstimateHess') && isa(o2, 'EstimateHess')
                newest = cell(size(o1.estimates,1));
                for ii = 1:size(o1.estimates,1)
                    for ij = 1:size(o1.estimates,2)
                        if isempty(o1.estimates{ii,ij})
                            newest{ii,ij} = o2.estimates{ii,ij};
                        elseif isempty(o2.estimates{ii})
                            newest{ii,ij} = o1.estimates{ii,ij};
                        else
                            newest{ii,ij} = o1.estimates{ii,ij} - o2.estimates{ii,ij};
                        end
                    end
                end
            elseif isa(o1, 'EstimateHess')
                newest = cell(size(o1.estimates,1));
                for ii = 1:size(o1.estimates,1)
                    for ij = 1:size(o1.estimates,2)
                        newest{ii,ij} = o1.estimates{ii,ij} - o2;
                    end
                end
            else
                newest = cell(size(o2.estimates,1));
                for ii = 1:size(o2.estimates,1)
                    for ij = 1:size(o2.estimates,2)
                        newest{ii,ij} = o1 - o2.estimates{ii,ij};
                    end
                end
            end
            r = EstimateHess(newest);
        end

        function r = uminus(o1)
            newest = cell(size(o1.estimates,1));
            for ii = 1:size(o1.estimates,1)
                for ij = 1:size(o1.estimates,2)
                    newest{ii,ij} = -o1.estimates{ii,ij};
                end
            end
            r = EstimateHess(newest);
        end

        function r = mtimes(o1, o2)
            if ~isa(o1, 'EstimateHess')
                error('Type error!');
            end
            if ~isa(o2, 'EstimateVec')
                error('Type error!')
            end
            Hess = o1.Matrix();
            grad = o2.Vector();
            tr = Hess*grad;
            r = EstimateVec(tr, o1.GetSizes());
        end

        function r = mldivide(o1, o2)
            if ~isa(o1, 'EstimateHess')
                error('Type error!');
            end
            Hess = o1.Matrix();
            if isa(o2, 'EstimateVec')
                grad = o2.Vector();
            else
                grad = o2;
            end
            tr = Hess\grad;
            r = EstimateVec(tr, o1.GetSizes());
        end

        function self = set.estimates(self, est)
            if size(est,1) ~= size(est,2)
                error('"estimates" property has to be a square cell array.');
            end
            for ii = 1:size(est,1)
                for ij = (ii+1):size(est,2)
                    if ~isempty(est{ij,ii})
                        if isempty(est{ii,ij})
                            est{ii,ij} = est{ij,ii}';
                        end
                        est{ij,ii} = [];
                        warning('All derivatives are stored in upper triangle part of the "estimates" cell array. Values in lower triangle part are in some cases lost!');
                    end
                end
            end
            self.estimates = est;
        end

        function sizes = GetSizes(self)
            nest = size(self.estimates,1);
            sizes = zeros(nest,1);
            for ii = 1:nest
                sizes(ii) = size(self.estimates{ii,ii},1);
            end
        end

        function H = Matrix(self)
            nest = size(self.estimates,1);
            sizes = self.GetSizes();
            H = zeros(sum(sizes));
            libound = 1;
            ljbound = 1;
            for ii = 1:nest
                H(libound:libound+sizes(ii)-1, ljbound:ljbound+sizes(ii)-1) = self.estimates{ii,ii};
                ljbound = ljbound + sizes(ii);
                for ij = ii+1:nest
                    H(libound:libound+sizes(ii)-1, ljbound:ljbound+sizes(ij)-1) = self.estimates{ii,ij};
                    H(ljbound:ljbound+sizes(ij)-1, libound:libound+sizes(ii)-1) = self.estimates{ii,ij}';
                    ljbound = ljbound + sizes(ij);
                end
                libound = libound + sizes(ii);
                ljbound = libound;
            end
        end

    end

end