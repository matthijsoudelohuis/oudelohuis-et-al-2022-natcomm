
params.Experiments          = {'VisOnlyTwolevels' 'ChangeDetectionConflict'};
params.area                 = 'V1';

%% Get data:
[Data]              = MOL_GetData('E:','CHDET',params.Experiments,{},[],{'sessionData' 'trialData' 'spikeData'});
sessionData         = Data.sessionData;
trialData           = Data.trialData;
spikeData           = Data.spikeData;

%% Remove last 20 trials:
trialData           = MOL_RemoveLastnTrials(trialData,20);

%% Filter out neurons based on quality:
spikeData           = MOL_filterNeurons(sessionData,trialData,spikeData);

%% Filter sessions with photostimulation in V1 and virus in V1:
idx                                 = sessionData.UseOpto & strcmp(sessionData.PhotostimArea,'V1') & strcmp(sessionData.OpsinArea,'V1');
[sessionData,trialData,spikeData]   = MOL_getTempPerSes(sessionData.session_ID(idx),sessionData,trialData,spikeData);

%% Filter neurons on area:
spikeFields     = fieldnames(spikeData);
idx             = ismember(spikeData.area,params.area);
for iF = 1:length(spikeFields)
    spikeData.(spikeFields{iF}) = spikeData.(spikeFields{iF})(idx);
end

%% Get only sessions with neurons recorded:
[sessionData,trialData,spikeData] = MOL_getTempPerSes(unique(spikeData.session_ID),sessionData,trialData,spikeData);
fprintf('Dataset: %d sessions, %d trials, %d neurons\n',length(sessionData.session_ID),length(trialData.session_ID),length(spikeData.session_ID));

%% Save the dataset:
save('E:\Matlab\MOL_Analysis\AOudeLohuisetal_2022_NatComms\7Noisecorrelations\Data7_1.mat')