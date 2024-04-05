classdef Estimate_vec
%This is a class used for estimating multiple parameters simultaneously,
%e.g. conductivity + contact impedances. estimates-property contains a cell
%array of different estimated parameters, e.g. estimates{1} = conductivity,
%estimates{2} = contact impedances.

    properties
        estimates
    end


    methods

        function obj = Estimate_vec(estimates, sizes)
            if nargin == 1
                obj.estimates = estimates;
            elseif nargin == 2
                est = cell(length(sizes),1);
                libound = 1;
                for ii = 1:length(sizes)
                    est{ii} = estimates(libound:libound+sizes(ii)-1);
                    libound = libound + sizes(ii);
                end
                obj.estimates = est;
            end
        end

        function r = plus(o1, o2)
            if isa(o1, 'Estimate_vec') && isa(o2, 'Estimate_vec')
                newest = cell(length(o1.estimates),1);
                for ii = 1:length(o1.estimates)
                    if isempty(o1.estimates{ii})
                        newest{ii} = o2.estimates{ii};
                    elseif isempty(o2.estimates{ii})
                        newest{ii} = o1.estimates{ii};
                    else
                        newest{ii} = o1.estimates{ii} + o2.estimates{ii};
                    end
                end
            elseif isa(o1, 'Estimate_vec')
                newest = cell(length(o1.estimates),1);
                for ii = 1:length(o1.estimates)
                    newest{ii} = o1.estimates{ii} + o2;
                end
            else
                newest = cell(length(o2.estimates),1);
                for ii = 1:length(o2.estimates)
                    newest{ii} = o2.estimates{ii} + o1;
                end
            end
            r = Estimate_vec(newest);
        end

        function r = minus(o1, o2)
            if isa(o1, 'Estimate_vec') && isa(o2, 'Estimate_vec')
                newest = cell(length(o1.estimates),1);
                for ii = 1:length(o1.estimates)
                    if isempty(o1.estimates{ii})
                        newest{ii} = o2.estimates{ii};
                    elseif isempty(o2.estimates{ii})
                        newest{ii} = o1.estimates{ii};
                    else
                        newest{ii} = o1.estimates{ii} - o2.estimates{ii};
                    end
                end
            elseif isa(o1, 'Estimate_vec')
                newest = cell(length(o1.estimates),1);
                for ii = 1:length(o1.estimates)
                    newest{ii} = o1.estimates{ii} - o2;
                end
            else
                newest = cell(length(o2.estimates),1);
                for ii = 1:length(o2.estimates)
                    newest{ii} = o1 - o2.estimates{ii};
                end
            end
            r = Estimate_vec(newest);
        end

        function r = uminus(o1)
            newest = cell(length(o1.estimates),1);
            for ii = 1:length(o1.estimates)
                newest{ii} = -o1.estimates{ii};
            end
            r = Estimate_vec(newest);
        end

        function r = mtimes(o1, o2)
            if isa(o1, 'Estimate_vec')
                newest = cell(length(o1.estimates),1);
                for ii = 1:length(o1.estimates)
                    newest{ii} = o1.estimates{ii}*o2;
                end
            elseif isa(o2, 'Estimate_vec')
                newest = cell(length(o2.estimates),1);
                for ii = 1:length(o2.estimates)
                    newest{ii} = o1*o2.estimates{ii};
                end
            end
            r = Estimate_vec(newest);
        end

        function r = norm(o1)
            r = 0;
            for ii = 1:length(o1.estimates)
                r = r + sum(o1.estimates{ii}.^2);
            end
            r = sqrt(r);
        end

        function vec = Vector(self)
            nest = size(self.estimates,1);
            sizes = zeros(nest,1);
            for ii = 1:nest
                sizes(ii) = size(self.estimates{ii},1);
            end
            vec = zeros(sum(sizes),1);
            libound = 1;
            for ii = 1:nest
                vec(libound:libound+sizes(ii)-1) = self.estimates{ii};
                libound = libound + sizes(ii);
            end
        end

        function sizes = GetSizes(self)
            nest = size(self.estimates,1);
            sizes = zeros(nest,1);
            for ii = 1:nest
                sizes(ii) = size(self.estimates{ii},1);
            end
        end



    end

end