function [A, C ,outparams] = gas_wrapper(data,varargin)
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

%%%% STARTING MESSAGES PART FOR THIS RUN
global VERBOSE LOGIT
VERBOSE = true;
LOGIT = true;

if nargin == 3
    params = varargin{1};
    gastype = varargin{2};
    arq_connect = struct();
else
    arq_connect = varargin{1};
    gastype = [];
end


if isfield(arq_connect, 'params')&&~isempty(arq_connect.params)
    params = arq_connect.params;
end
if isfield(arq_connect, 'method')&&~isempty(arq_connect.method)
    gastype = arq_connect.method;
elseif 0
    gastype = 'gwr';
    %error('Method must be declared')
end
if ~isfield(params, 'savegas')
    params.savegas.name = 'gas';
    params.savegas.resume = false;
    params.savegas.path = '~/Dropbox/octave_progs/';
end

if isempty(params)
    
    NODES = 100;
    
    params = struct();
    
    params.use_gpu = false; %%% the use of gpu is very very very poor, so it runs much slower than it should. do not use
    params.PLOTIT = true; %
    params.RANDOMSTART = false; % if true it overrides the .startingpoint variable
    params.RANDOMSET = false;
    params.savegas.name = 'gas';
    params.savegas.resume = false;
    params.savegas.path = '~/Dropbox/octave_progs/';
    
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
else
    %
    
end
if  isfield(arq_connect, 'name')&&params.savegas.resume
    params.savegas.name = strcat(arq_connect.name,'-n', num2str(params.nodes), '-s',num2str(size(data,1)),'-q',num2str(params.q),'-i',num2str(labindex));
elseif params.savegas.resume
    error('Strange arq_connect definition. ''.name'' field is needed.')
end


MAX_EPOCHS = params.MAX_EPOCHS;
PLOTIT = params.PLOTIT;

%%% things that are specific for skeletons:
if isfield(params,'skelldef')
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
if ~isfield(params,'use_gpu')|| gpuDeviceCount==0
    params.use_gpu = false;
    %or break or error or warning...
end

if PLOTIT
    figure
    plotgwr() % clears plot variables
end

datasetsize = size(data,2);
errorvect = nan(1,MAX_EPOCHS*datasetsize);
epochvect = nan(1,MAX_EPOCHS*datasetsize);
nodesvect = nan(1,MAX_EPOCHS*datasetsize);

if  params.use_gpu
    data = gpuArray(data);
    errorvect = gpuArray(errorvect);
    epochvect = gpuArray(epochvect);
    nodesvect = gpuArray(nodesvect);
end

if ~isfield(params, 'savegas')||~isfield(params.savegas, 'parallelgasescount')||(isfield(params.savegas, 'parallelgases')&&~params.savegas.parallelgases)||~isfield(params.savegas, 'parallelgases')
    params.savegas.parallelgasescount = 0;
    dbgmsg('Warning: no parallelgases variable. update caller code!',1)
end

gasindex = params.savegas.parallelgasescount+labindex;

if isfield(params, 'savegas')&&isfield(params.savegas, 'parallelgasescount')&&isfield(params.savegas, 'parallelgases')&&params.savegas.parallelgases
    gasmaxcount = params.savegas.P+params.savegas.parallelgasescount;
else
    gasmaxcount = 1;
end

[gasgas,gasfun] = create_gas(gastype,params,data);


%%% ok, I will need to be able to resume working on a gas if I want to
if isfield(params, 'savegas')&&params.savegas.resume
    savegas = strcat(params.savegas.path,'/', params.savegas.name,'.mat');
    %%% checks if savegas exists
    disp('checkpoint1')
    if exist(savegas,'file') %%% if it does, then it loads it, because it should have the same
        dbgmsg('Found gas with the same characteristics as this one. Will try loading gas',params.savegas.name,1)
        load(savegas)
        disp('checkpoint2')
        try
            gasfun(data(:,1), gasgas(gasindex,end));
            
            gasacuu = zeros(1,size(gasgas,2));
            for iii = 1:size(gasgas,2)
                gasacuu(1,iii) = gasgas(gasindex,iii).params.accumulatedepochs;
            end
            if isfield(params.savegas, 'accurate_track_epochs')&& params.savegas.accurate_track_epochs
                [epochsuntilnow,bestgasindex] = max(gasacuu(gasacuu<MAX_EPOCHS));
            else
                [epochsuntilnow,bestgasindex] = max(gasacuu);
            end
            disp('checkpoint3')
            if ~isempty(bestgasindex)
                gasgas(gasmaxcount,end+1) = gas;
                gasgas(gasindex,end) = gasgas(gasindex,bestgasindex);
                MAX_EPOCHS = MAX_EPOCHS-epochsuntilnow;
            else
                if size(gasgas,1)<gasmaxcount
                    gasgas(gasmaxcount,end+1) = gas;
                    gasgas(gasindex,end) = create_gas(gastype,params,data);
                else
                    gasgas(gasindex,end+1) = create_gas(gastype,params,data);
                end
                disp('checkpoint4')
            end
            %                 if gasgas(iii).params.accumulatedepochs > MAX_EPOCHS
            %                     dbgmsg('saved gas exceeds the amount of epochs, have to start new gas');
            %                     error('saved gas exceeds the amount of epochs, have to start new gas') %oh, I shouldnt be doing this
            %                 end
            %   
            
        catch
            dbgmsg('I failed. Will start a new gas. ',1)
            if isfield(params.savegas, 'accurate_track_epochs')&& ~params.savegas.accurate_track_epochs&&~isfield(gasgas(gasindex,1).params, 'accumulatedepochs')
                %%% this should load old gases done before
                %%% accurate_track_epochs was created.
                gasgas(1:gasmaxcount,end+1) =gasgas(gasindex,1);
            elseif gasindex==1
                gasgas(1:gasmaxcount,end+1) = create_gas(gastype,params,data);
            else
                gasgas(gasindex,end) = create_gas(gastype,params,data);
            end
        end
    end
    
else
    dbgmsg('Didn''t find gas with the same characteristics as this one. Will use a new gas with name:',params.savegas.name,1)
end

%%%%% enabling saved parallel gases is a pain
if isfield(params, 'savegas')&&params.savegas.resume



therealk = 0; %% a real counter for epochs

%%%starting main loop

for num_of_epochs = 1:MAX_EPOCHS % strange idea: go through the dataset more times - actually this makes it overfit the data, but, still it is interesting.
    
    if params.RANDOMSET
        kset = randperm(datasetsize);
    else
        kset = 1:datasetsize;
    end
    % start of the loop
    disp(gasindex)
    disp(size(gasgas))
    
    for k = kset %step 1
        therealk = therealk +1;
        
        gasgas(gasindex,end) = gasfun(data(:,k), gasgas(gasindex,end));
        
        %to make it look nice...
        errorvect(therealk) = gasgas(gasindex,end).a;
        epochvect(therealk) = therealk;
        nodesvect(therealk) = gasgas(gasindex,end).r;
        if PLOTIT&&mod(k,plottingstep)==0&&numlabs==1 %%% also checks to see if it is inside a parpool
            plotgwr(gasgas(gasindex,end).A,gasgas(gasindex,end).C,errorvect,epochvect,nodesvect, skelldef, layertype)
            drawnow
        end
    end
    
    %%% now save the resulting gas
    if isfield(params, 'savegas')&&params.savegas.resume
        save(savegas, 'gasgas')
    end
end

%%% updating the number of epochs the gas has run
gasgas(gasindex,end) = gasgas(gasindex,end).update_epochs(MAX_EPOCHS);

%%% now save the end-resulting gas
if isfield(params, 'savegas')&&params.savegas.resume
    dbgmsg('Saving gas:',params.savegas.name,1)
    save(savegas, 'gasgas')
end
outparams.graph.errorvect = errorvect;
outparams.graph.epochvect = epochvect;
outparams.graph.nodesvect = nodesvect;
outparams.accumulatedepochs = gasgas(gasindex,end).params.accumulatedepochs;
outparams.initialnodes = [gasgas(gasindex,end).ni1,gasgas(gasindex,end).ni2];
A = gasgas(gasindex,end).A;
C = gasgas(gasindex,end).C;
end
end
function [gasgas,gasfun] = create_gas(gastype,params,data)
switch gastype
    case 'gwr'
        gasfun = @gwr_core;
        gasgas = gas;
        gasgas = gasgas.gwr_create(params,data);
    case 'gng'
        gasfun = @gng_core;
        gasgas = gas;
        gasgas = gasgas.gng_create(params,data);
    otherwise
        error('Unknown method.')
end
end