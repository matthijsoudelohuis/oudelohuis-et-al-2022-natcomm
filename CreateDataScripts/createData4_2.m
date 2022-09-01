%% Reset all
startover

%% Parameter settings:
params                      = params_histresponse_coding; %Parameters for PSTH (All time is in microseconds)
params.conv_win             = 'chg';

params.Experiments          = {'ChangeDetectionConflictDecor' 'VisOnlyTwolevels' 'ChangeDetectionConflict' }; %Which versions of the task to load data from
params.ExperimentLabels     = {'NE' 'UST' 'MST'}; %Labels for each task version
params.nExperiments         = length(params.Experiments);

params.AlignOn              = 'stimChange';      %On which timestamp to align as t=0

params.area                 = 'V1'; %Filter only V1 data

params                      = MOL_getColors_CHDET(params);

params.savedir              = 'E:\Documents\PhD\Figures\Project CHDET\Results - 2nd Bump\2 Neural\Zmat_AlignStim';

params.nshuffle             = 1000;
params.alpha                = 0.05;

params.minNneurons          = 10;
params.markersize           = 100;

params.minTrialCond         = 10;

%% Get data:
[Data] = MOL_GetData('E:','CHDET',params.Experiments,{},[],{'sessionData' 'trialData_newtrials' 'spikeData'});
sessionData     = Data.sessionData;
trialData       = Data.trialData;
spikeData       = Data.spikeData;
sessionData.Experiment = strrep(sessionData.Experiment,'ChangeDetectionConflict2','ChangeDetectionConflict');

%% Remove last 20 trials:
trialData = MOL_RemoveLastnTrials(trialData,20);

%% Filter neurons on area:
spikeFields     = fieldnames(spikeData);
idx             = ismember(spikeData.area,params.area);
for iF = 1:length(spikeFields)
    spikeData.(spikeFields{iF}) = spikeData.(spikeFields{iF})(idx);
end
fprintf('Filtered neurons based on area\n');

%% Filter out multisensory sessions that have too low visual or auditory performance:
nSessions           = length(sessionData.session_ID);
visperf             = NaN(nSessions,1);
auperf              = NaN(nSessions,1);
for iSes = 1:nSessions
    sesidx          = strcmp(trialData.session_ID,sessionData.session_ID(iSes));
    vistrialidx     = strcmp(trialData.trialType,'X') & trialData.visualOriChangeNorm==3 & trialData.hasphotostim~=1 & sesidx; 
    visperf(iSes)   = sum(trialData.correctResponse(vistrialidx)) / sum(vistrialidx); 
    autrialidx      = strcmp(trialData.trialType,'Y') & trialData.audioFreqChangeNorm==3 & trialData.hasphotostim~=1 & sesidx;
    auperf(iSes)    = sum(trialData.correctResponse(autrialidx)) / sum(autrialidx); 
end

sesids              = sessionData.session_ID(~(strcmp(sessionData.Experiment,'ChangeDetectionConflict') & (visperf<0.3 | auperf<0.3)));
fprintf('Removed %d/%d sessions with low behavioral accuracy\n',length(sessionData.session_ID)-length(sesids),length(sessionData.session_ID));
[sessionData, trialData, spikeData]        = MOL_getTempPerSes(sesids,sessionData,trialData,spikeData);

%% Filter out sessions with muscimol:
sesids              = sessionData.session_ID(cellfun(@isempty,sessionData.MuscimolArea));
fprintf('Removed %d/%d sessions  with muscimol manipulations\n',length(sessionData.session_ID)-length(sesids),length(sessionData.session_ID));
[sessionData, trialData, spikeData]        = MOL_getTempPerSes(sesids,sessionData,trialData,spikeData);

%% Filter out trials with photostimulation in V1:
% sesids              = sessionData.session_ID(~(sessionData.UseOpto & (strcmp(sessionData.PhotostimArea,'V1') | strcmp(sessionData.PhotostimArea,'PPC'))));
sesids              = sessionData.session_ID(sessionData.UseOpto & strcmp(sessionData.PhotostimArea,'V1'));
idx                 = trialData.hasphotostim==1 & ismember(trialData.session_ID,sesids);
trialFields         = fieldnames(trialData);
fprintf('Removed %d/%d trials  with photostim in V1\n',sum(idx),length(trialData.session_ID));
for iF = 1:length(trialFields)
    trialData.(trialFields{iF}) = trialData.(trialFields{iF})(~idx);
end

%% Filter out trials with reaction times <200ms;
idx                 = trialData.responseLatency<200e3;
trialFields         = fieldnames(trialData);
for iF = 1:length(trialFields)
    trialData.(trialFields{iF}) = trialData.(trialFields{iF})(~idx);
end

%% Get only sessions with neurons recorded:
[sessionData,trialData,spikeData] = MOL_getTempPerSes(unique(spikeData.session_ID),sessionData,trialData,spikeData);
fprintf('Dataset: %d sessions, %d trials, %d neurons\n',length(sessionData.session_ID),length(trialData.session_ID),length(spikeData.session_ID));

%%




%% Set parameters and init output fields:

params.zscore               = 1;
params.smoothing            = 1;
params.binsize              = 10e3;             %Size of the bins
params.conv_win             = 'chg';       %Type of window used for smoothing {flat, gaussian)
params.conv_sigma           = 50e3;           %sd of gaussian window for smoothing

%Construct bin edges and time axis
params.edges                = [params.t_pre:params.binsize:params.t_post] - params.binsize/2;                    %#ok<NBRAK> Define edges of the bins
params.xtime                = params.t_pre:params.binsize:params.t_post-params.binsize;                        %Define time axis
params.nTimebins            = length(params.xtime); %number of time bins

nContrasts                  = 2;

nNeurons                    = length(spikeData.ts);
lastsesid                   = []; %memory var to keep track of whether neuron comes from different session and design matrix needs to be reconstructed
fprintf('Computing average firing rate for neuron        \n');

%Some measures of dimensionality:
uSessions                   = unique(sessionData.session_ID);
nSessions                   = length(uSessions);

outputmat                   = NaN(nSessions,params.nTimebins,nContrasts); %init output variables:
outputmat_shuf              = NaN(nSessions,params.nTimebins,nContrasts,params.nshuffle);
outputmat_n                 = NaN(nSessions,1);

%% For each session compute AUC
fprintf('Computing responses for session     \n');


for iSes = 1:nSessions
    fprintf('Computing AUC for session %d/%d\n',iSes,nSessions)
    sesid       = sessionData.session_ID(iSes);
    %Get the relevant data for each session individually:
    [temptrialData,tempspikeData] = MOL_getTempPerSes(sesid,trialData,spikeData);
    
    %Construct tensor:
    events_ts       = temptrialData.(params.AlignOn);
    nNeurons        = length(tempspikeData.ts);
    nTimebins       = length(params.xtime);
    nTrials         = length(events_ts);
    tensor          = NaN(nNeurons,nTimebins,nTrials);
    
    for iNeuron = 1:nNeurons
        %Get histogram per neuron
        hist_mat                = calc_psth(events_ts,tempspikeData.ts{iNeuron},params);
        tensor(iNeuron,:,:)     = hist_mat';
    end
    
    splits                  = {};
    splits{1}               = ismember(temptrialData.trialType,{'X'}) & temptrialData.visualOriChangeNorm==2 & temptrialData.vecResponse==3;
    splits{2}               = ismember(temptrialData.trialType,{'X'}) & temptrialData.visualOriChangeNorm==2 & temptrialData.vecResponse==2;
    splits{3}               = ismember(temptrialData.trialType,{'X'}) & temptrialData.visualOriChangeNorm==3 & temptrialData.vecResponse==3;
    splits{4}               = ismember(temptrialData.trialType,{'X'}) & temptrialData.visualOriChangeNorm==3 & temptrialData.vecResponse==2;
    
    popmat                  = squeeze(nanmean(tensor,1));
    
    contrasts               = [1 2; 3 4];
    
    outputmat_n(iSes,1)     = nNeurons;
    
%     Compute AUC:
    for iContr = 1:nContrasts %for each contrast (split vs split trial sets)
        if sum(splits{contrasts(iContr,1)})>params.minTrialCond && sum(splits{contrasts(iContr,2)})>params.minTrialCond
            for iTime = 1:params.nTimebins %loop over all timebins
                responses                               = [popmat(iTime,splits{contrasts(iContr,1)})'; popmat(iTime,splits{contrasts(iContr,2)})']; %get all firing rates in this bin for all trials
                feature                                 = [ones(sum(splits{contrasts(iContr,1)}),1); ones(sum(splits{contrasts(iContr,2)}),1)*2]; %Make correct feature vector (ones for condition A and two for condition B)
                [~,~,~,outputmat(iSes,iTime,iContr)]    = perfcurve(feature,responses,1);
                %Make shuffle distribution and compute permutation test threshold:
                AUC_shuffle = NaN(1,params.nshuffle); %init vector for storing AUC shuffled values
                for it = 1:params.nshuffle
                    [~,~,~,AUC_shuffle(it)]    = perfcurve(feature(randperm(length(feature))),responses,1);
                end
                outputmat_shuf(iSes,iTime,iContr,:) = AUC_shuffle;
            end
        end
    end
end

%% Save the data:
save('E:\Matlab\MOL_Analysis\AOudeLohuisetal_2022_NatComms\4AUC\Data4_2.mat','spikeData','trialData','sessionData','outputmat','outputmat','outputmat_shuf','outputmat_n')



