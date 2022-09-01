%% Verification_EarlyLateSilencing
% This script aligns V1 neural activity to stimulus onset and compare averaged 
% firing rate of V1 neurons during optogenetic silencing with early or late onset (0 versus 200 ms)
% Reproduces figure 5d
% Oude Lohuis et al. 2022 Nat Comms

%% Parameter settings for PSTH
params                      = params_histresponse(); % All time is in microseconds
params.zscore               = 0; %no z-scoring
params.conv_sigma           = 0.01e6;        %sd of gaussian window for smoothing %smaller than for the rest of the analysis
params.conv_twin            = 9*params.conv_sigma;         %Window size for smoothing

%% SnakePlot Parameters:
params.SortBy               = 'maxResponse';
params.AlignOn              = 'stimChange';      %On which timestamp to align as t=0

params.clipzscore           = 20;               %remove neurons that have  
params.clipmodrate          = 1.5;              %remove neurons that are 
params.normBaseline         = 1;                %Whether firing rate is normalized to baseline (100%)
params.cscale               = [-100 300];

params.area                 = 'V1';

%% Load the data:
load('E:\Matlab\MOL_Analysis\AOudeLohuisetal_2022_NatComms\5Opto\Data5_2.mat')
fprintf('Dataset: %d sessions, %d trials, %d neurons\n',length(sessionData.session_ID),length(trialData.session_ID),length(spikeData.session_ID));

%% Base splits:
params.nSplits                  = 3;
params.colors_splits            = {[0 0 0.6] [0.7 0.7 1] [0 0 0.8]};
params.labels_splits            = {'Control' 'Early Inactivation' 'Late inactivation'};

%% Main loop to get psth matrix:
nNeurons                = length(spikeData.ts);
lastsesid               = []; %memory var to keep track of whether neuron comes from different session and design matrix needs to be reconstructed
fprintf('Computing firing rate for neuron        \n');
snakemat                = NaN(nNeurons,params.nTimebins,params.nSplits);

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
    splits{1}       = strcmp(temptrialData.trialType,'X') & temptrialData.visualOriChangeNorm==3 & ...
        temptrialData.hasphotostim==0 & temptrialData.vecResponse==2;
    splits{2}       = strcmp(temptrialData.trialType,'X') & temptrialData.visualOriChangeNorm==3 &...
        temptrialData.hasphotostim==1 & temptrialData.PostChangeOptoStart==0 & (temptrialData.optoEnd-temptrialData.stimChange)>0.8e6;
    splits{3}       = strcmp(temptrialData.trialType,'X') & temptrialData.visualOriChangeNorm==3 & ...
        temptrialData.hasphotostim==1 & temptrialData.PostChangeOptoStart==0.2 & (temptrialData.optoEnd-temptrialData.stimChange)>0.8e6;
    
    for iSplit = 1:params.nSplits
        hist_mean                       = nanmean(hist_mat(splits{iSplit},:),1);
        %Normalize:
        if params.normBaseline
            hist_mean                   = hist_mean/nanmean(hist_mean(params.xtime<0));
        end
        snakemat(iNeuron,:,iSplit)      = hist_mean;
    end
end

%% Sort snakemat
if isfield(params,'SortBy')
    params.twin_resp_start  = 0.05e6;
    params.twin_resp_stop   = 0.2e6;
    %     [maxresp,maxidx]          = max(snakemat_splits{1}(:,params.xtime>params.twin_resp_start & params.xtime<=params.twin_resp_stop),[],2);
    [maxresp,maxidx]          = max(snakemat(:,params.xtime>params.twin_resp_start & params.xtime<=params.twin_resp_stop,1),[],2);
    
    if strcmp(params.SortBy,'peakLatency')
        [~,sortidx]             = sort(maxidx,1,'descend');
        snakemat                = snakemat(sortidx,:,:);
    elseif strcmp(params.SortBy,'maxResponse')
        [~,sortidx]             = sort(maxresp,1,'descend');
        snakemat                = snakemat(sortidx,:,:);
    end
end

%% Remove clipped firing rates (impossible values):
idx_z         = any(any(snakemat>params.clipzscore,3),2);
if any(idx_z)
    snakemat    = snakemat(~idx_z,:,:);
    fprintf('Removed %d neurons with unreasonable z-scores\n',sum(idx_z));
end

%% Remove activated neurons:
modulationratio = mean(snakemat(:,params.xtime>0 & params.xtime<1e6,2),2) ./ mean(snakemat(:,params.xtime>-1e6 & params.xtime<0,2),2);
idx_m = modulationratio>params.clipmodrate;
if any(idx_m)
    snakemat    = snakemat(~idx_m,:,:);
    fprintf('Removed %d neurons because of increase in modulation rate\n',sum(idx_m))
end

%% Make figure of the zscored matrix in an imagesc:
figure; set(gcf,'units','normalized','Position',[0.1 0.1 0.6 0.5],'color','w')
for iSplit = 1:length(splits)
    subplot(1,length(splits),iSplit)
    imagesc(snakemat(:,:,iSplit)*100,params.cscale); hold on;
    plot([find(params.xtime == 0) find(params.xtime == 0)], [0 size(snakemat,1)+0.5],'k','LineWidth',5);
    set(gca, 'XTick', 1:1000:length(params.xtime), 'XTickLabels', params.xtime(1:1000:length(params.xtime))/1e6,'FontSize', 20)
    set(gca, 'YTick', [1 size(snakemat,1)], 'YTickLabels', [1 size(snakemat,1)],'FontSize', 20)
    xlabel('Time (s)','FontSize', 20)
    ylabel('Neuron','FontSize', 20)
    c = colorbar;
    colormap(parula);
    c.Label.String = 'Firing rate (% of baseline)';
end

%% Make figure of the mean over time:
figure; set(gcf,'units','normalized','Position',[0.3 0.2 0.6 0.5],'color','w')
handles = [];
for iSplit = 1:length(splits)
    meantoplot = nanmean(snakemat(:,:,iSplit),1) * 100;
    errortoplot = nanstd(snakemat(:,:,iSplit),1)/sqrt(size(snakemat,1)) * 100;
    h = shadedErrorBar(params.xtime,meantoplot,errortoplot,{'-k','markerfacecolor',params.colors_splits{iSplit},'LineWidth',3},1);
    h.mainLine.Color = params.colors_splits{iSplit};    h.patch.FaceColor = params.colors_splits{iSplit};
    h.edge(1).Color = [1 1 1];  h.edge(2).Color =  [1 1 1];
    h.edge(1).LineWidth = 0.1;  h.edge(2).LineWidth = 0.1;
    handles(iSplit) = h.mainLine; hold all; %#ok<SAGROW>
end

set(gca, 'XTick', [-0.3e6 0 0.2e6 0.5e6 1e6 1.5e6], 'XTickLabels', [-0.3e6 0 0.2e6 0.5e6 1e6 1.5e6]./1e6,'FontSize', 20)
xlim([-0.5e6 1.5e6]);

ylabel('Firing rate (norm to baseline)','FontSize', 20)
xlabel('Time (s)','FontSize', 20)
legend(handles,params.labels_splits); legend boxoff

