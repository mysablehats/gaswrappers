function [A, C ,outparams] = gas_wrapper(data,params,gastype)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%cf parisi, 2015 and cf marsland, 2002
%based on the GNG algorithm from the guy that did the GNG algorithm for
%matlab

% some tiny differences:
% in the 2002 paper, they want to show the learning of topologies ability
% of the GWR algorithm, which is not our main goal. In this sense they have
% a function that can generate new points as pleased p(eta). This is not
% our case, we will just go through our data sequentially

% I am not taking time into account. the h(time) function is therefore
% something that yields a constant value

%the initial parameters for the algorithm:
%global maxnodes at en eb h0 ab an tb tn amax
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%maxnodes = params.nodes; %maximum number of nodes/neurons in the gas
%at = params.at;%0.95; %activity threshold
%en = params.en;%= 0.006; %epsilon subscript n
%eb = params.eb;%= 0.2; %epsilon subscript b
%amax = params.amax;%= 50; %greatest allowed age
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

if isempty(params)
    
    NODES = 100;
    
    params = struct();
    
    params.PLOTIT = true; %not really working
    params.RANDOMSTART = false; % if true it overrides the .startingpoint variable
    
    n = randperm(size(data,2),2);
    params.startingpoint = [n(1) n(2)];
    
    params.amax = 500; %greatest allowed age
    params.nodes = NODES; %maximum number of nodes/neurons in the gas
    params.en = 0.006; %epsilon subscript n
    params.eb = 0.2; %epsilon subscript b
    
    %Exclusive for gwr
    params.STATIC = true;
    params.MAX_EPOCHS = 1; % this means data will be run over twice
    params.at = 0.80; %activity threshold
    params.h0 = 1;
    params.ab = 0.95;
    params.an = 0.95;
    params.tb = 3.33;
    params.tn = 3.33;
    
    %Exclusive for gng
    params.age_inc                  = 1;
    params.lambda                   = 3;
    params.alpha                    = .5;     % q and f units error reduction constant.
    params.d                           = .99;   % Error reduction factor.

end

MAX_EPOCHS = params.MAX_EPOCHS;
PLOTIT = params.PLOTIT;

%%% things that are specific for skeletons:
if isfield(params,'skeldef')
    skelldef = params.skelldef;
else
    skelldef = [];
end
if isfield(params,'layertype')
    layertype = params.layertype;
else
    layertype = [];
end

if isfield(params,'plottingstep')        
    if params.plottingstep == 0
        plottingstep = size(data,2);
    else
        plottingstep = params.plottingstep;
    end    
else
    plottingstep = fix(size(data,2)/20);
end

if PLOTIT
    figure
    plotgwr() % clears plot variables
end

datasetsize = size(data,2);
errorvect = nan(1,MAX_EPOCHS*datasetsize);
epochvect = nan(1,MAX_EPOCHS*datasetsize);
nodesvect = nan(1,MAX_EPOCHS*datasetsize);

switch gastype
    case 'gwr'
        gasfun = @gwr_core;
        gasgas = gas;
        gasgas = gasgas.gwr_create(params,data);
    case 'gng'
        gasfun = @gng_core;
        gasgas = gas;
        gasgas = gasgas.gng_create(params,data);
end

therealk = 0; %% a real counter for epochs

for num_of_epochs = 1:MAX_EPOCHS % strange idea: go through the dataset more times - actually this makes it overfit the data, but, still it is interesting.
    
    % start of the loop
    for k = 1:datasetsize %step 1
        therealk = therealk +1;
        
        gasgas = gasfun(data(:,k), gasgas);
        
        %to make it look nice...
        errorvect(therealk) = gasgas.a;
        epochvect(therealk) = therealk;
        nodesvect(therealk) = gasgas.r;
        if PLOTIT&&mod(k,plottingstep)==0
            plotgwr(gasgas.A,gasgas.C,errorvect,epochvect,nodesvect, skelldef, layertype)
            drawnow
        end
        
    end
end
outparams.graph.errorvect = errorvect;
outparams.graph.epochvect = epochvect;
outparams.graph.nodesvect = nodesvect;
outparams.initialnodes = [gasgas.ni1,gasgas.ni2];
A = gasgas.A;
C = gasgas.C;
end
