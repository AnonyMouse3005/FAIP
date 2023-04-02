% Facility location w/ penalty object for LP-based algorithm from Charikar 2001
% Remark: because penalties in our characterization are not large enough,
% the algorithms are more willing to deny services to clients.
% Hence, some clients are not connected to any facility at any time
classdef FLP < handle
    properties  % class properties
        z (1,1) double {mustBeNonnegative};
        t (1,1) double {mustBeNonnegative};  % current time
        cost (:,:) double {mustBeNonnegative};
        penalty (1,:) double {mustBeNonnegative};
        alpha (1,:) double {mustBeNonnegative};
    end
    properties (SetAccess = protected)  % all these attributes change in response to alpha
        % client nodes
        connwit (1,:) double {mustBeNonnegative};  % connecting witness for each client
        connected (1,:) double {mustBeNumericOrLogical};  % whether alpha_j is connected i.e., still growable
        timeout (1,:) double {mustBeNumericOrLogical};
        % facility nodes
        t_paid (1,:) double;  % time at which facility i gets paid for (Inf initially)
        temp_open (1,:) double {mustBeNumericOrLogical};  % whether facility i is temporarily open
        % edges
        beta (:,:) double {mustBeNonnegative};  % dual variables for each edge between facility i and client j
        tight (:,:) double {mustBeNumericOrLogical};  % whether edge (i,j) is tight
    end
    methods
        function obj = FLP(v,n,z,cost,penalty)  % constructor
            obj.z = z;
            obj.t = 1;
            obj.cost = cost;
            obj.penalty = penalty;
            obj.connwit = zeros(1,n);
            obj.connected = false(1,n);
            obj.timeout = false(1,n);
            obj.t_paid = inf(1,v);
            obj.temp_open = false(1,v);
            obj.tight = false(v,n);
            obj.beta = zeros(v,n);
            obj.alpha = zeros(1,n);
        end
        function set.alpha(obj, value)  % update beta whenever alpha changes
            obj.alpha = value;
            obj.freezeAlpha();
            obj.declareTight();
        end
        function set.tight(obj, value)  % connect j (if unconnected) to i whenever edge (i,j) is tight and i is temp. open
            obj.tight = value;
            obj.growBeta_conn2Fac();
        end
        function set.beta(obj, value)  % update facilities' status whenever beta changes
            obj.beta = value;
            obj.openFacility();
        end
        function growAlpha(obj)
            J = find(~obj.connected & ~obj.timeout);
            % choose alpha from an unconnected client (LEXICOGRAPHICALLY) to grow
            for j = 1:length(J)
                obj.alpha(J(j)) = obj.alpha(J(j)) + 1;  % grow alpha uniformly
            end
        end
    end
    methods (Access = private)
        function declareTight(obj)
            J = find(~obj.connected);
            for j = 1:length(J)
                % set of closed cities i in which edge (i,j) is tight
                I = find(obj.alpha(J(j)) >= obj.cost(:,J(j)) & ~obj.temp_open');
                if ~isempty(I)  % pick one arbitrary city to declare (i,j) tight
%                     obj.tight(randsample(I,1),J(j)) = true;
                    obj.tight(I(1),J(j)) = true;
                end
            end
        end
        function growBeta_conn2Fac(obj)
            J = find(~obj.connected);
            for j = 1:length(J)
                I = find(obj.tight(:,J(j)));
                for i = 1:length(I)  % edge (i,j) is tight
                    if ~obj.temp_open(I(i))  % consider growing beta(i,j) if j is not connected and i is not open
                        % beta(i,j) + cost(i,j) must not exceed alpha(j)
                        if obj.beta(I(i),J(j)) + obj.cost(I(i),J(j)) < obj.alpha(J(j))
                            obj.beta(I(i),J(j)) = obj.beta(I(i),J(j)) + 1;
                        end
                    else  % connect j to i if j is not connected but i is already open
                        obj.connected(J(j)) = true;
                        obj.connwit(J(j)) = i;
                    end
                end
            end
        end
        function openFacility(obj)
            I = find(sum(obj.beta,2) >= obj.z & ~obj.temp_open');  % all closed facility i's that are paid for
            for i = 1:length(I)
                [~, J] = find(obj.tight(I(i),:) & ~obj.connected);  % (1 x ?) all unconnected client j's that have a tight edge to i
                J = unique(J);
                for j = 1:length(J)
                    obj.temp_open(I(i)) = true;  % temporarily open i
                    obj.t_paid(I(i)) = obj.t;
                    obj.connected(J(j)) = true;
                    obj.connwit(J(j)) = I(i);
                end
            end
        end
        function freezeAlpha(obj)
            J = find(obj.alpha == obj.penalty & ~obj.connected);
            obj.timeout(J) = true;
        end
    end
end
