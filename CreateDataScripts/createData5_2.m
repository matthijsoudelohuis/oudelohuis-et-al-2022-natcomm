params.area                 = 'V1';

%% Get data:
[Data] = MOL_GetData('E:','CHDET',{'VisOnlyTwolevels' 'ChangeDetectionConflict'},{},[],{'sessionData' 'trialData' 'spikeData'});
sessionData     = Data.sessionData;
trialData       = Data.trialData;
spikeData       = Data.spikeData;

%% Filter out neurons based on quality:
spikeData       = MOL_filterNeurons(sessionData,trialData,spikeData);

%% Filter sessions:
idx                                 = sessionData.UseOpto & strcmp(sessionData.PhotostimArea,'V1') & strcmp(sessionData.OpsinArea,'V1');
fprintf('Selected %d/%d sessions optogenetic stimulation in V1 with recordings in V1\n',sum(idx),length(sessionData.session_ID));
[sessionData,trialData,spikeData]   = MOL_getTempPerSes(sessionData.session_ID(idx),sessionData,trialData,spikeData);

%% Filter neurons on area:
spikeFields     = fieldnames(spikeData);
idx             = ismember(spikeData.area,params.area);

for iF = 1:length(spikeFields)
    spikeData.(spikeFields{iF}) = spikeData.(spikeFields{iF})(idx);
end

%% Filter out trials with reaction times <200ms; 
idx                 = trialData.responseLatency<200e3;
trialFields         = fieldnames(trialData);
for iF = 1:length(trialFields)
    trialData.(trialFields{iF}) = trialData.(trialFields{iF})(~idx);
end

%% save the data:
save('E:\Matlab\MOL_Analysis\AOudeLohuisetal_2022_NatComms\5Opto\Data5_2.mat','sessionData','trialData','spikeData')
