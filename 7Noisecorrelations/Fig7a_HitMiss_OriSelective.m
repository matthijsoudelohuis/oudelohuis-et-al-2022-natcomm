%% This script shows the averaged firing rate for hits and misses for all orientation selective neurons in V1
% Reproduces fig 7a
% Oude Lohuis et al. 2022 Nat Comms

%% Reset all
startover

%% Load dataset:
load('E:\Matlab\MOL_Analysis\AOudeLohuisetal_2022_NatComms\4AUC\Data4_1.mat')
fprintf('Dataset: %d sessions, %d trials, %d neurons\n',length(sessionData.session_ID),length(trialData.session_ID),length(spikeData.session_ID));

%% Parameter settings:
params.Experiments          = {'ChangeDetectionConflictDecor' 'VisOnlyTwolevels' 'ChangeDetectionConflict' }; %Which versions of the task to load data from
params.ExperimentLabels     = {'NE' 'UST' 'MST'}; %Labels for each task version
params.nExperiments         = length(params.Experiments);

params.area                 = 'V1'; %Filter only V1 data

params                      = MOL_getColors_CHDET(params);

params.minTrialCond         = 3;

% parameters for auROC:
params.nBinThreshold        = 5;

%% Main loop to get psth matrix:
params.nSplits          = 4;
nNeurons                = length(spikeData.session_ID);
lastsesid               = []; %memory var to keep track of whether neuron comes from different session and design matrix needs to be reconstructed
fprintf('Computing average Z-scored response for neuron        \n');
zmat                    = NaN(nNeurons,params.nTimebins,params.nSplits);

for iNeuron = 1:nNeurons %Loop over all neurons:
    fprintf(repmat('\b', 1, numel([num2str(iNeuron-1) num2str(nNeurons)])+2));
    fprintf('%d/%d\n',iNeuron,nNeurons);
    
    if ~strcmp(lastsesid,spikeData.session_ID(iNeuron)) %construct new predictor matrix if neuron comes from a new session:
        %Get the relevant data for each session individually:
        temptrialData        = MOL_getTempPerSes(spikeData.session_ID(iNeuron),trialData);
        lastsesid            = spikeData.session_ID(iNeuron); %save this session_ID
    end
    
    %Compute histogram:
    events_ts               = temptrialData.(params.AlignOn);
    hist_mat                = calc_psth(events_ts,spikeData.ts{iNeuron},params);    %Construct histogram matrix:
    
    splits          = {};
    splits{1}       = ismember(temptrialData.trialType,{'X'}) & temptrialData.visualOriChangeNorm==3 & ismember(temptrialData.visualOriPostChangeNorm,[1 2]) & temptrialData.vecResponse==3 & ~(temptrialData.hasphotostim==1);
    splits{2}       = ismember(temptrialData.trialType,{'X'}) & temptrialData.visualOriChangeNorm==3 & ismember(temptrialData.visualOriPostChangeNorm,[1 2]) & temptrialData.vecResponse==2 & ~(temptrialData.hasphotostim==1);
    splits{3}       = ismember(temptrialData.trialType,{'X'}) & temptrialData.visualOriChangeNorm==3 & ismember(temptrialData.visualOriPostChangeNorm,[3 4]) & temptrialData.vecResponse==3 & ~(temptrialData.hasphotostim==1);
    splits{4}       = ismember(temptrialData.trialType,{'X'}) & temptrialData.visualOriChangeNorm==3 & ismember(temptrialData.visualOriPostChangeNorm,[3 4]) & temptrialData.vecResponse==2 & ~(temptrialData.hasphotostim==1);
    
    for iSplit = 1:params.nSplits
        if sum(splits{iSplit})>=params.minTrialCond
            zmat(iNeuron,:,iSplit) = mean(hist_mat(splits{iSplit},:),1);
        end
    end
end

%% Show the average firing rate

%Take timebins during first 200ms:
timeidx             = params.xtime>0 & params.xtime<0.2e6;
%Identify neurons that have at least one significant bin during this first 200ms:
idx_signori         = sum(outputmat_sign(2,:,timeidx),3)>=1;
%Take all the neurons from UST and MST:
iExp                = [2 3];
idx_exp             = ismember(spikeData.session_ID,sessionData.session_ID(ismember(sessionData.Experiment,params.Experiments(iExp))));

%Reorder firing rate matrix to preferred orientation of selective units:
zmatoriresp                     = squeeze(nanmean(zmat(:,timeidx,:),2));
pref_ori_12                     = zmatoriresp(:,1) >= zmatoriresp(:,3);

zmat2                           = zmat;
zmat2(~pref_ori_12,:,[1 2])     = zmat(~pref_ori_12,:,[3 4]);
zmat2(~pref_ori_12,:,[3 4])     = zmat(~pref_ori_12,:,[1 2]);


idx                             = idx_signori' & idx_exp; %Combine indices
datatoplot                      = squeeze(nanmean(zmat2(idx,:,:),1)); %average over neurons
errortoplot                     = squeeze(nanstd(zmat2(idx,:,:),1)); %std over neurons

params.splitcolors              = {[0.1 0.3 0.1] [0.1 0.3 0.1] [0.1 0.1 0.3] [0.1 0.1 0.3]};
params.splitlines               = {'-.' '-' '-.' '-'};

%Make the figure:
figure; set(gcf,'units','normalized','Position',[0.2 0.1 0.4 0.55],'color','w'); hold all;
title('Orientation coding','FontSize',10)

handles = NaN(params.nSplits,1);
for iSplit = 1:params.nSplits
    meantoplot = squeeze(nanmean(zmat2(idx,:,iSplit),1));
    errortoplot = squeeze(nanstd(zmat2(idx,:,iSplit),1)) / sqrt(sum(idx));
    h = shadedErrorBar(params.xtime,meantoplot,errortoplot,{params.splitlines{iSplit},'markerfacecolor',params.splitcolors{iSplit},'LineWidth',4},1);
    h.mainLine.Color = params.splitcolors{iSplit};    h.patch.FaceColor = params.splitcolors{iSplit};
    delete(h.edge(1)); delete(h.edge(2));
    handles(iSplit) = h.mainLine; hold all;
end

%Figure make-up:
plot([0 0],[0 8],'k:','LineWidth',2)
set(gca, 'XTick', [-0.5e6 0 0.5e6 1e6], 'XTickLabels', [-0.5e6 0 0.5e6 1e6]/1e6,'FontSize', 20)
xlim([-0.2e6 1e6]);
set(gca, 'YTick', [4 12] , 'FontSize', 20,'FontName','Arial')
ylim([2 12]);
ylabel('Firing rate (sp/s)','FontSize', 20)
xlabel('Time (s)','FontSize', 20)
legend(handles,{'Pref Ori - Miss' 'Pref Ori - Hit' 'Nonpref Ori - Miss' 'Nonpref Ori - Hit'}); legend boxoff


%%