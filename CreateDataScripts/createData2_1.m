params.Experiments          = {'ChangeDetectionConflictDecor' 'VisOnlyTwolevels' 'ChangeDetectionConflict' }; %Which versions of the task to load data from
params.area                 = 'V1'; %Filter only V1 data

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

%% Filter out neurons based on quality:
spikeData       = MOL_filterNeurons(sessionData,trialData,spikeData);

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

%% Save the dataset:
save('E:\Matlab\MOL_Analysis\AOudeLohuisetal_2022_NatComms\2Neural\Data2_1.mat','sessionData','trialData','spikeData')
