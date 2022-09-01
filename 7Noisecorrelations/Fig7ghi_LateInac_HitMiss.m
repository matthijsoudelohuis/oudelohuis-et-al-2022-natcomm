%% This script analyzes the V1 data just preceding late inactivation
% It analyzes the firing rate and the noise correlations in the 0-200 ms window
% following visual stimuli preceding a hit or a miss
% Oude Lohuis et al. 2022 Nat Comms

%% Parameter settings for PSTH
params                      = params_histresponse_coding(); % All time is in microseconds
params.zscore               = 1;
params.smoothing            = 1;
params.conv_sigma           = 0.01e6;        %sd of gaussian window for smoothing
params.conv_twin            = round((9*params.conv_sigma)/params.binsize)*params.binsize;         %Window size for smoothing

%E.g. method can be either large squared bins or 1 ms bins with smoothing window ('smoothing') or combination thereof
params.conv_win             = 'gaussian';

% parameters for window of interest:
params.t_pre                = -1e6;
params.t_post               = 0.2e6;

%Construct bin edges and time axis
params.edges                = [params.t_pre:params.binsize:params.t_post] - params.binsize/2;                    %#ok<NBRAK> Define edges of the bins
params.xtime                = params.t_pre:params.binsize:params.t_post-params.binsize;                        %Define time axis
params.nTimebins            = length(params.xtime); %number of time bins

params.AlignOn              = 'stimChange';      %On which timestamp to align as t=0

params.area                 = 'V1';

params.Experiments          = {'VisOnlyTwolevels' 'ChangeDetectionConflict'};
params.ExperimentLabels     = {'UST' 'MST'};
params.nExperiments         = length(params.Experiments);

params.nSplits              = 4;

params                      = MOL_getColors_CHDET(params);

%Statistical level of signficiance: 
params.alpha                = 0.05;

%% Load the dataset:
load('E:\Matlab\MOL_Analysis\AOudeLohuisetal_2022_NatComms\7Noisecorrelations\Data7_1.mat')

%% Main loop to compute firing rate:
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
    splits{1}       = strcmp(temptrialData.trialType,'X') & temptrialData.visualOriChangeNorm==3 ...
        & temptrialData.hasphotostim==1 & temptrialData.PostChangeOptoStart==0.2 & temptrialData.vecResponse==2;
    splits{2}       = strcmp(temptrialData.trialType,'X') & temptrialData.visualOriChangeNorm==3 ...
        & temptrialData.hasphotostim==1 & temptrialData.PostChangeOptoStart==0.2 & temptrialData.vecResponse==3;
    splits{3}       = strcmp(temptrialData.trialType,'X') & temptrialData.visualOriChangeNorm==2 ...
        & temptrialData.hasphotostim==1 & temptrialData.PostChangeOptoStart==0.2 & temptrialData.vecResponse==2;
    splits{4}       = strcmp(temptrialData.trialType,'X') & temptrialData.visualOriChangeNorm==2 ...
        & temptrialData.hasphotostim==1 & temptrialData.PostChangeOptoStart==0.2 & temptrialData.vecResponse==3;

    for iSplit = 1:params.nSplits
        snakemat(iNeuron,:,iSplit) = nanmean(hist_mat(splits{iSplit},:),1);
    end
end

% %% Make figure of the mean:
% figure; set(gcf,'units','normalized','Position',[0.3 0.2 0.26 0.4],'color','w')
% 
% idx_all = true(size(spikeData.session_ID));
% 
% handles = [];
% for iSplit = 1:params.nSplits
%     meantoplot      = nanmean(snakemat(idx_all,:,iSplit),1);
%     errortoplot     = nanstd(snakemat(idx_all,:,iSplit),1)/sqrt(sum(idx_all));
% 
%     h = shadedErrorBar(params.xtime,meantoplot,errortoplot,{params.lines_splits{iSplit},'markerfacecolor',params.colors_splits{iSplit},'LineWidth',3},1);
%     h.mainLine.Color = params.colors_splits{iSplit};    h.patch.FaceColor = params.colors_splits{iSplit};
%     delete(h.edge(1)); delete(h.edge(2)); %h.edge(1).Color = [1 1 1];  h.edge(2).Color =  [1 1 1];
% %     h.edge(1).LineWidth = 0.1;  h.edge(2).LineWidth = 0.1;
%     handles(iSplit) = h.mainLine; hold all; %#ok<SAGROW>
% end
% 
% set(gca, 'XTick', [-0.2e6 0 0.2e6 0.4e6 1e6 1.5e6], 'XTickLabels', [-0.2e6 0 0.2e6 0.4e6 1e6 1.5e6]./1e6,'FontSize', 20)
% xlim([-0.2e6 0.4e6]);
% 
% ylabel('Z-scored firing rate','FontSize', 20)
% xlabel('Time (s)','FontSize', 20)
% legend(handles,params.labels_splits); legend boxoff


%% Make figure of the mean:

params.lines_splits             = {'-' '--' '-' '--'};
params.labels_splits            = {'Max - Hit' 'Max - Miss' 'Thr - Hit' 'Thr - Miss'};

params.colors_splits       = {  params.colors_experiments{2} params.colors_experiments{2} params.colors_experiments{2}.^[0.4 0.2 1] params.colors_experiments{2}.^[0.4 0.2 1];
                                params.colors_experiments{3} params.colors_experiments{3} params.colors_experiments{3}.^[0.4 1 1] params.colors_experiments{3}.^[0.4 1 1]};

figure; set(gcf,'units','normalized','Position',[0.3 0.2 0.36 0.24],'color','w'); hold all;
for iExp = 1:2
    subplot(1,2,iExp); hold all;
    idx_exp = ismember(spikeData.session_ID,sessionData.session_ID(ismember(sessionData.Experiment,params.Experiments{iExp})));
    handles = [];

    for iSplit = 1:params.nSplits
        meantoplot      = nanmean(snakemat(idx_exp,:,iSplit),1);
        handles(iSplit) = plot(params.xtime,meantoplot,params.lines_splits{iSplit},'Color',params.colors_splits{iExp,iSplit},'LineWidth',3); %#ok<SAGROW>
    end
    %figure makeup:
    plot([0 0],[-0.2 0.3],':','Color',[0.6 0.6 0.6],'LineWidth',2)
    patch([200e3 250e3 250e3 200e3],[-0.2 -0.2 0.6 0.6],[0 0.6 0.85],'EdgeColor',[0 0.6 0.85],'LineWidth',3)
    set(gca, 'XTick', [-0.2e6 0 0.2e6 0.4e6 1e6 1.5e6], 'XTickLabels', [-0.2e6 0 0.2e6 0.4e6 1e6 1.5e6]./1e6,'FontSize', 10)
    xlim([-0.2e6 0.25e6]);
    ylim([-0.2 0.6])
    ylabel('Z-scored firing rate','FontSize', 12)
    xlabel('Time (s)','FontSize', 12)
    idx = [1 0 1 0];
    legend(handles(idx==1),params.labels_splits(idx==1),'Location','NorthWest'); legend boxoff
end

%% Make bar figure of the hit-miss difference in 100-200ms: 

snakematdiff                = []; %Construct difference between hits and misses:
snakematdiff(:,:,1)         = diff(snakemat(:,:,[2 1]),[],3); %Diff Max Hit and Miss
snakematdiff(:,:,2)         = diff(snakemat(:,:,[4 3]),[],3);  %Diff Thr Hit and Miss

params.colors_splits       = {  params.colors_experiments{2} params.colors_experiments{2}.^[0.4 0.2 1];
                                params.colors_experiments{3} params.colors_experiments{3}.^[0.4 1 1]};

datatoplot = NaN(2,2,200);
idx_time = params.xtime>100e3 & params.xtime<190e3; %Do not include last 10ms for jitter effects of photostimulation on resp

for iExp = 1:2
    idx_exp         = ismember(spikeData.session_ID,sessionData.session_ID(ismember(sessionData.Experiment,params.Experiments{iExp})));
    for iSplit = 1:2
        datatoplot(iExp,iSplit,1:sum(idx_exp))      = squeeze(nanmean(snakematdiff(idx_exp,idx_time,iSplit),2));
    end
end

meantoplot          = nanmean(datatoplot,3);
errortoplot         = nanstd(datatoplot,[],3)/sqrt(sum(idx_exp));

figure; set(gcf,'units','normalized','Position',[0.3 0.2 0.3 0.4],'color','w'); hold all;
handles = [];
for iExp = 1:2
    for iSplit = 1:2
        barloc = (iSplit-1)*2+iExp;
        handles(end+1)   = bar(barloc,meantoplot(iExp,iSplit),0.8);
        set(handles(end),'facecolor',params.colors_splits{iExp,iSplit});
        errorbar(barloc,meantoplot(iExp,iSplit),errortoplot(iExp,iSplit),'k.')
    end
end

%% Statistical testing:

respmat         = squeeze(nanmean(snakemat(:,idx_time,:),2));

G_mou           = cell(nNeurons,1);
uMice           = unique(sessionData.mousename);
for iMouse = 1:length(uMice)
    G_mou(ismember(spikeData.session_ID,sessionData.session_ID(strcmp(sessionData.mousename,uMice{iMouse})))) = uMice(iMouse);
end

for iExp = 1:2
    idx_exp         = ismember(spikeData.session_ID,sessionData.session_ID(ismember(sessionData.Experiment,params.Experiments{iExp})));
    fprintf('%s:      \n',params.ExperimentLabels{iExp})
    for iSplit = 1:2
        Y_resp = [respmat(:,1+(iSplit-1)*2); respmat(:,2+(iSplit-1)*2)];
        X_hitm = [ones(nNeurons,1); ones(nNeurons,1)*2];
        
        tbl             = table(Y_resp([idx_exp; idx_exp]),X_hitm([idx_exp; idx_exp]),[G_mou(idx_exp); G_mou(idx_exp)],'VariableNames',{'Activity','HitMiss','Mouse'}); %Create table for mixed model
        lme             = fitlme(tbl,'Activity~HitMiss+(1|Mouse)'); %construct linear mixed effects model with fixed effect of temporal window and random intercept for different mice
        stats           = dataset2table(anova(lme,'DFMethod','Satterthwaite')); %Perform ANOVA on model and output as matrix
        fprintf('%s vs Miss: ',params.labels_splits{iSplit + (iSplit-1)})
        fprintf('F(%d,%2.0f)=%1.2f, p=%1.3f; \n',stats{2,3},stats{2,4},stats{2,2},stats{2,5})
        p = stats{2,5};
        if p<0.05
            barloc = (iSplit-1)*2+iExp;
            sigstar([barloc-0.1 barloc+0.1],p)
        end
        
        tempfile = fullfile(sprintf('SourceData_Fig7g_Activity_HitMiss_%s_%s.xlsx',params.ExperimentLabels{iExp}, params.labels_splits{iSplit*2-1}(1:3)));
        writetable(tbl,tempfile)
    end
end

%%






%% Compute noise correlations:
params.binsize              = 25e3;
params.conv_sigma           = 50e3;
params.conv_twin            = round((9*params.conv_sigma)/params.binsize)*params.binsize;         %Window size for smoothing
params.conv_win             = 'gaussian';
params.minNtrialscond       = 10;

%Construct bin edges and time axis
params.edges                = [params.t_pre:params.binsize:params.t_post] - params.binsize/2;                    %#ok<NBRAK> Define edges of the bins
params.xtime                = params.t_pre:params.binsize:params.t_post-params.binsize;                        %Define time axis
params.nTimebins            = length(params.xtime); %number of time bins

params.nSplits = 2;

%Some measures of dimensionality:
nTotalNeurons               = length(spikeData.session_ID);
nSessions                   = length(sessionData.session_ID);

noisecorrmat                = NaN(params.nSplits,nTotalNeurons,nTotalNeurons,params.nTimebins); %init output matrix

neuroncounter               = 0;

for iSes = 1:nSessions
    fprintf('Computing noise correlations for session %d/%d\n',iSes,nSessions)
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
    
    splits          = {};
    splits{1}       = strcmp(temptrialData.trialType,'X') ...
        & temptrialData.hasphotostim==1 & temptrialData.PostChangeOptoStart==0.2 & temptrialData.vecResponse==2;
    splits{2}       = strcmp(temptrialData.trialType,'X') ...
        & temptrialData.hasphotostim==1 & temptrialData.PostChangeOptoStart==0.2 & temptrialData.vecResponse==3;
    
    %subtract mean response for each condition:
    for iSplit = 1:params.nSplits
        tensor(:,:,splits{iSplit})    = tensor(:,:,splits{iSplit}) - repmat(nanmean(tensor(:,:,splits{iSplit}),3),1,1,sum(splits{iSplit}));
    end
    
    %Compute noise corr:
    for iSplit = 1:params.nSplits

        if sum(splits{iSplit})>params.minNtrialscond
            
            for iNx = 1:nNeurons
                for iNy = 1:nNeurons
                    for iT = 1:params.nTimebins
                        noisecorrmat(iSplit,neuroncounter + iNx,neuroncounter + iNy,iT) = corr(squeeze(tensor(iNx,iT,splits{iSplit})),squeeze(tensor(iNy,iT,splits{iSplit})));
                    end
                end
            end
            
        end
    end
    
    neuroncounter = neuroncounter + length(tempspikeData.session_ID);
end

%Set all correlation values of 1 to NaN (autocorrelation)
noisecorrmat(noisecorrmat>0.99)   = NaN;

%% Get all noise corr values according to condition and reshape:
nCross              = nTotalNeurons^2;
datatoplot          = NaN(params.nSplits,nCross,params.nTimebins);

for iSplit = 1:params.nSplits
    temp                        = noisecorrmat(iSplit,:,:,:);
    temp                        = reshape(temp,nCross,params.nTimebins);
    datatoplot(iSplit,:,:)      = temp;
end

idx             = all(all(~isnan(datatoplot),3),1);
fprintf('n = %d pairwise combinations\n',sum(idx))

%Compute mean and sem:
meantoplot      = squeeze(nanmean(datatoplot,2));
errortoplot     = squeeze(nanstd(datatoplot,[],2) ./ sqrt(sum(~isnan(datatoplot),2)));

%% Figure of noise correlations over time for hit and miss:
figure; set(gcf,'units','normalized','Position',[0.3 0.2 0.26 0.4],'color','w'); hold all;
handles = [];

for iSplit = 1:params.nSplits
    h = shadedErrorBar(params.xtime, meantoplot(iSplit,:),squeeze(errortoplot(iSplit,:,:)),{params.lines_splits{1,iSplit},'LineWidth',3,'Color',params.colors_splits{iSplit}},0);
    delete(h.edge(1)); delete(h.edge(2)); 
    handles(end+1) = h.mainLine; %#ok<SAGROW>
end

statstoplot     = false(params.nTimebins,1);

idx             = all(all(~isnan(datatoplot),3),1);

for iT = 1:params.nTimebins
    m = bootstrp(1000,@mean,diff(datatoplot(:,idx,iT),[],1));
    if prctile(m,1.25)>0 || prctile(m,98.75)<0
        statstoplot(iT) = true;
    end
end

for iT = 2:params.nTimebins-1
    if statstoplot(iT)
        plot([params.xtime(iT) params.xtime(iT+1)],[0.115 0.115],'Color','k','LineWidth',5)
    end
end

plot([0 0],[0 0.2],'k','LineWidth',2)
ylabel('Average noise correlation','FontSize',20)
legend(handles,{'Hit' 'Miss'}); legend boxoff;
xlim([-0.2e6 0.4e6]);
ylim([0.04 0.12])
xlabel('Time (ms)')
set(gca,'XTick',-0.2e6:0.2e6:1e6,'XTickLabels',(-0.2e6:0.2e6:1e6)*1e-3,'YTick',[0.04 0.08 0.12])

tblHit = array2table(squeeze(datatoplot(1,:,:)),'VariableNames',sprintfc('Timepoint_%g',1:48));
tblMiss = array2table(squeeze(datatoplot(2,:,:)),'VariableNames',sprintfc('Timepoint_%g',1:48));

tempfile = fullfile('SourceData_Fig7i_NCoverTime_HitMiss.xlsx');
writetable(tblHit,tempfile,'Sheet',1)
writetable(tblMiss,tempfile,'Sheet',2)

%% Jean: removed lines with NaN for lisibility in excel

tblHit2 = table2array(tblHit);
tblMiss2 = table2array(tblMiss);
tblHit3 = tblHit2(~isnan(nanmean(tblHit2,2)), :);
tblMiss3 = tblMiss2(~isnan(nanmean(tblMiss2,2)), :);
writematrix(tblHit3,'SourceData_Fig7i_JP.xlsx','Sheet',1)
writematrix(tblMiss3,'SourceData_Fig7i_JP.xlsx','Sheet',2)




