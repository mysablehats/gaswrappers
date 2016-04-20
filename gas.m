classdef gas
    properties
        params
        A
        C
        C_age
        hizero
        hszero
        time
        t0
        r
        h
        a
        awk
        ni1
        ni2
    end
    methods
        function X = hs(gasgas, t)
            h0 = gasgas.params.h0;%= 1;
            ab = gasgas.params.ab;%= 0.95;
            tb = gasgas.params.tb;%= 3.33;
            X = h0 - gasgas.S(t)/ab*(1-exp(-ab*t/tb));
        end
        function X = hi(gasgas, t)
            h0 = gasgas.params.h0;%= 1;
            an = gasgas.params.an;%= 0.95;
            tn = gasgas.params.tn;%= 3.33;
            X = h0 - gasgas.S(t)/an*(1-exp(-an*t/tn));
        end
        function X = S(~,t)
            X = 1;
        end
        function [n1,n2, ni1,ni2] = initialnodes(gasgas, data)
            
            RANDOMSTART = gasgas.params.RANDOMSTART;
            
            % (1)
            % pick n1 and n2 from data
            ni1 = 1;
            ni2 = 2;
            if RANDOMSTART
                n = randperm(size(data,2),2);
                ni1 = n(1);
                ni2 = n(2);
            elseif isfield(gasgas.params,'startingpoint')&&~isempty('params.startingpoint')
                ni1 = gasgas.params.startingpoint(1);
                ni2 = gasgas.params.startingpoint(2);
            end
            if ni1> size(data,2)
                dbgmsg('n1 overflow. using last data point',num2str(ni1),1)
                n1 = data(:,end);
            else
                n1 = data(:,ni1);
            end
            if ni2> size(data,2)
                dbgmsg('n2 overflow. using last data point',num2str(ni2),1)
                n2 = data(:,end);
            else
                n2 = data(:,ni2);
            end
            
        end
        function gasgas = gas_create(gasgas, params,data)
            gasgas.params = params;
            
            [n1,n2, gasgas.ni1,gasgas.ni2] = initialnodes(gasgas, data);
            
            gasgas.A = zeros(size(n1,1),params.nodes);
            gasgas.A(:,[1 2]) = [n1, n2];
            
            %%%%%% NEW ADDITION: the adjusting weighting constant parameter
            %%%% this is not so simple because the input vectors
            if isfield(gasgas.params, 'skelldef')&&isfield(gasgas.params.skelldef, 'awk')&&isfield(params, 'layertype')
                %%% setting awk based on which layer I am in
                %%% inputlayers, q(1) and layertype influence awk :'(
                switch params.layertype
                    case 'pos'
                        gasgas.awk = repmat(params.skelldef.awk.pos,1,params.q(1));
                    case 'vel'
                        gasgas.awk = repmat(params.skelldef.awk.vel,1,params.q(1));
                end
                try
                    n1.*gasgas.awk;
                catch
                    dbgmsg('Tried to use your awk definition, but I failed.',1)
                    gasgas.awk = ones(size(n1,1),1);
                end
            else
                gasgas.awk = ones(size(n1,1),1);
            end
        end
        function gasgas = gwr_create(gasgas, params, data)
            
            gasgas = gasgas.gas_create(params,data);
            
            %%%gwr specific initialization
            %test some algorithm conditions:
            if ~(0 < params.en || params.en < params.eb || params.eb < 1)
                error('en and/or eb definitions are wrong. They are as: 0<en<eb<1.')
            end
            
            % (2)
            % initialize empty set C
            
            gasgas.C = sparse(params.nodes,params.nodes); % this is the connection matrix.
            gasgas.C_age = gasgas.C;
            
            gasgas.r = 3; %the first point to be added is the point 3 because we already have n1 and n2
            gasgas.h = zeros(1,params.nodes);%firing counter matrix
            
            %%% SPEEDUP CHANGE
            if params.STATIC
                gasgas.hizero = gasgas.hi(0)*ones(1,params.nodes);
                gasgas.hszero = gasgas.hs(0);
            else
                gasgas.time = 0;
            end
            gasgas.t0 = cputime; % my algorithm is not necessarily static!
        end
    end
end