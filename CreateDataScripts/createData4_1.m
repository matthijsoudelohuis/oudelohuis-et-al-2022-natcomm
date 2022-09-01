%% First load the auROC results:
loaddate                = '09-08-20';
datadir                 = fullfile('E:\Data\Analysis\auROC\',loaddate);
nVars                   = 33;

load(fullfile(datadir,'auROC_1.mat'))
outputmat               = NaN(nVars,size(outputmat,1),size(outputmat,2)); %#ok<*NODEF>
outputmat_thr           = NaN(nVars,size(outputmat_thr,1),size(outputmat_thr,2));
outputmat_nSpikes       = NaN(nVars,size(outputmat_nSpikes,1),size(outputmat_nSpikes,2));

for iVar = 1:nVars
    loadstruct                      = load(fullfile(datadir,sprintf('auROC_%d.mat',iVar)),'outputmat','outputmat_thr','outputmat_nSpikes');
    %     loadstruct                      = load(fullfile(datadir,sprintf('auROC_%d.mat',iVar)),'outputmat','outputmat_thr');
    outputmat(iVar,:,:)             = loadstruct.outputmat;
    outputmat_thr(iVar,:,:,:)       = loadstruct.outputmat_thr;
    outputmat_nSpikes(iVar,:,:,:)   = loadstruct.outputmat_nSpikes;
end

params.xtime_AUC            = params.xtime;
params.area                 = 'V1'; %Filter only V1 data

%% Parameter that defines which bins have a sufficient number of spikes to base permutation test on.
params.nBinThreshold        = 5;
outputmat(outputmat_nSpikes<params.nBinThreshold) = NaN;

%% compute significance;
outputmat_sign      = outputmat > outputmat_thr;

%% Filter neurons on area:
spikeFields     = fieldnames(spikeData);
idx             = ismember(spikeData.area,params.area);
for iF = 1:length(spikeFields)
    spikeData.(spikeFields{iF}) = spikeData.(spikeFields{iF})(idx);
end
outputmat           = outputmat(:,idx,:);
outputmat_sign      = outputmat_sign(:,idx,:);
outputmat_nSpikes   = outputmat_nSpikes(:,idx,:);
outputmat_thr       = outputmat_thr(:,idx,:);
fprintf('Selected only %d V1 neurons\n\n',sum(idx));

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
save('E:\Matlab\MOL_Analysis\AOudeLohuisetal_2022_NatComms\4AUC\Data4_1.mat')

