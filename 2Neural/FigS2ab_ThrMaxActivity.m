%% This script analyzes stimulus locked visually evoked activity for the 
% two different saliency of visual changes (threshold and maximal saliency)
% Supplementary Figure to:
% Oude Lohuis et al. 2022 Nat Comms

%% Reset all
startover

%% Parameter settings:
params                      = params_histresponse_coding; %Parameters for PSTH (All time is in microseconds)
params.conv_win             = 'chg';
params.zscore               = 1;

params.Experiments          = {'ChangeDetectionConflictDecor' 'VisOnlyTwolevels' 'ChangeDetectionConflict' }; %Which versions of the task to load data from
params.ExperimentLabels     = {'NE' 'UST' 'MST'}; %Labels for each task version
params.nExperiments         = length(params.Experiments);

params.AlignOn              = 'stimChange';      %On which timestamp to align as t=0

%parameters for display of zscore matrix:
params.colormap             = 'parula';
params.cscale               = [-0.5 2.2];

params.minTrialCond         = 3;

params.area                 = 'V1'; %Filter only V1 data

params.twin_resp_start      = 0e6; %only for sorting neurons
params.twin_resp_stop       = 0.5e6;

params                      = MOL_getColors_CHDET(params);

params.trialcategories      = 'ThrMax';
params.colors_ztrials       = {[156 111 240] [156 111 240] [2 27 148] [2 27 148]}; 
params.colors_ztrials       = cellfun(@(x) x/256,params.colors_ztrials,'UniformOutput',false);
params.lines_ztrials        = {'--' '-' '--' '-'};
params.labels_ztrials       = {'THR - MISS' 'THR - HIT' 'MAX - MISS' 'MAX - HIT'};

params.lines_experiments    = {'-' '-' '-'};

% params.savedir              = 'E:\Documents\PhD\Figures\Project CHDET\Results - 2nd Bump\2 Neural\Zmat_AlignStim';
params.savedir              = 'E:\Documents\PhD\Figures\Project CHDET\Results - 2nd Bump\2 Neural\Zmat_ThrMax';

%% Load dataset:
load('E:\Matlab\MOL_Analysis\AOudeLohuisetal_2022_NatComms\2Neural\Data2_1.mat')
fprintf('Dataset: %d sessions, %d trials, %d neurons\n',length(sessionData.session_ID),length(trialData.session_ID),length(spikeData.session_ID));

%% Show raster plot of example neurons showing visual responses to threshold and to maximal changes in orientation
cell_IDs            = {};
%Naive examples:
cell_IDs{end+1}     = '10092019030831134';
% cell_IDs{end+1}     = '10082019030831 218';
% cell_IDs{end+1}     = '10092019030721318';

% UST examples:
cell_IDs{end+1}     = '20292019121211218';

%MST examples:
% cell_IDs{end+1}     = '20032018020821278';
% cell_IDs{end+1}     = '20122018081431146';

% cell_IDs{end+1}     = '20122018081531198';
% cell_IDs{end+1}     = '20122018081531204';
% cell_IDs{end+1}     = '20222019062811204';
cell_IDs{end+1}     = '20122018081431146';
cell_IDs{end+1}     = '20092018082331113';

cell_IDs_zlabels        = {'a' 'b' 'c' 'd'};

outputmat_sign          = [];
params.AUC_varselec     = [];

MOL_plotRaster(sessionData,trialData,spikeData,cell_IDs,params)

%% Main loop to get psth matrix:
params.nSplits          = 4;

nNeurons                = length(spikeData.ts);
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
    
    splits                  = {};
    splits{1}               = ismember(temptrialData.trialType,{'X'}) & temptrialData.visualOriChangeNorm==2 & temptrialData.vecResponse==3;
    splits{2}               = ismember(temptrialData.trialType,{'X'}) & temptrialData.visualOriChangeNorm==2 & temptrialData.vecResponse==2;
    splits{3}               = ismember(temptrialData.trialType,{'X'}) & temptrialData.visualOriChangeNorm==3 & temptrialData.vecResponse==3;
    splits{4}               = ismember(temptrialData.trialType,{'X'}) & temptrialData.visualOriChangeNorm==3 & temptrialData.vecResponse==2;
    
    for iSplit = 1:params.nSplits %Store the mean response for each of these splits
        if sum(splits{iSplit})>=params.minTrialCond
            zmat(iNeuron,:,iSplit) = mean(hist_mat(splits{iSplit},:),1);
        end
    end
end

%% Sort zmat
meanresp                = mean(zmat(:,params.xtime>params.twin_resp_start & params.xtime<=params.twin_resp_stop,3),2) + mean(zmat(:,params.xtime>params.twin_resp_start & params.xtime<=params.twin_resp_stop,4),2);

[~,sortidx]             = sort(meanresp,1,'descend');

zmat                    = zmat(sortidx,:,:);            %Actually sort the zmats
spikeFields             = fieldnames(spikeData);    %Sort also the spikeData to match
for iF                  = 1:length(spikeFields)
    spikeData.(spikeFields{iF}) = spikeData.(spikeFields{iF})(sortidx,:);
end

%% Remove neurons with too few trials to compute averaged activity reliably:
idx_nan       = any(any(isnan(zmat),2),3);

%% Make figure:
for iExp = 1:params.nExperiments
    idx_exp = ismember(spikeData.session_ID,sessionData.session_ID(strcmp(sessionData.Experiment,params.Experiments(iExp))));
        
    idx_all = idx_exp & ~idx_nan;
    
    figure; set(gcf,'units','normalized','Position',[0.1 0.1 0.52 0.8],'color','w')
    suptitle(params.ExperimentLabels{iExp})
    for iSplit = 1:length(splits)
        subplot(2,2,iSplit)
        imagesc(zmat(idx_all,:,iSplit),params.cscale+0.1); hold on;
        plot([find(params.xtime == 0) find(params.xtime == 0)], [0 size(zmat(idx_all,:,iSplit),1)+0.5],'k:','LineWidth',3);
        
        set(gca, 'XTick', find(ismember(params.xtime,-3e6:0.5e6:3e6)), 'XTickLabels', params.xtime(ismember(params.xtime,(-3e6:0.5e6:3e6)))/1e6,'FontSize', 20)
        set(gca, 'YTick', [1 sum(idx_all)], 'YTickLabels', {1 sum(idx_all)},'FontSize', 20)
        xlim([find(params.xtime == -0.5e6) find(params.xtime == 1.5e6)]);
        xlabel('Time (s)','FontSize', 20)
        ylabel('Neuron','FontSize', 20)
        colormap(params.colormap);
        
        %Show example neurons on right y axis:
        yyaxis right
        set(gca, 'YTick',[])
        for i = 1:length(cell_IDs)
            temp = find(ismember(spikeData.cell_ID(idx_all),cell_IDs(i)));
            if temp
                yt = get(gca,'YTick'); yt(end+1) = sum(idx_all)-temp;
                ytl = get(gca,'YTickLabels'); ytl(end+1) = cell_IDs_zlabels{i};
                set(gca, 'YTick', yt, 'YTickLabels', ytl','FontSize', 20)
                ylim([1 sum(idx_all)])
            end
        end
    end
end

%make separate colorbar:
figure; set(gcf,'units','normalized','Position',[0.5 0.3 0.13 0.2],'color','w')
imagesc(zmat(idx_all,:,iSplit),params.cscale); hold on;
c = colorbar;
c.Label.String = 'Z-scored firing rate';
cla

%% For LMM statistics:
X_coh           = cell(nNeurons,1);
X_coh(ismember(spikeData.session_ID,sessionData.session_ID(ismember(sessionData.Experiment,params.Experiments{1}))))         = params.ExperimentLabels(1);
X_coh(ismember(spikeData.session_ID,sessionData.session_ID(ismember(sessionData.Experiment,params.Experiments{2}))))         = params.ExperimentLabels(2);
X_coh(ismember(spikeData.session_ID,sessionData.session_ID(ismember(sessionData.Experiment,params.Experiments{3}))))         = params.ExperimentLabels(3);

G_mou           = cell(nNeurons,1);
uMice           = unique(sessionData.mousename);
for iMouse = 1:length(uMice)
    G_mou(ismember(spikeData.session_ID,sessionData.session_ID(strcmp(sessionData.mousename,uMice{iMouse})))) = uMice(iMouse);
end

%% Make figure of the mean:
timeidx = params.xtime>0e3 & params.xtime<200e3;
fprintf('Activity difference per saliency: ')
for iExp = 1:params.nExperiments
    idx_exp = ismember(spikeData.session_ID,sessionData.session_ID(strcmp(sessionData.Experiment,params.Experiments(iExp))));
    idx_all = idx_exp;
    
   
    figure; set(gcf,'units','normalized','Position',[0.2 0.3 0.32 0.38],'color','w')
    suptitle(params.ExperimentLabels{iExp})
    
    %     handles = NaN(params.nSplits,1);
    handles = [];
    for iSplit = 1:params.nSplits
        meantoplot = nanmean(zmat(idx_all,:,iSplit),1);
        errortoplot = nanstd(zmat(idx_all,:,iSplit),1)/sqrt(size(zmat(idx_all,:,iSplit),1));
        plot([0 0], [-2 5],'k:','LineWidth',3);

        if ~all(isnan(meantoplot))
            h = shadedErrorBar(params.xtime,meantoplot,errortoplot,{params.lines_ztrials{iSplit},'markerfacecolor',params.colors_ztrials{iSplit},'LineWidth',3},0);
            h.mainLine.Color = params.colors_ztrials{iSplit};    h.patch.FaceColor = params.colors_ztrials{iSplit};
            delete(h.edge(1)); delete(h.edge(2)); %h.edge(1).Color = [1 1 1];  h.edge(2).Color =  [1 1 1]; h.edge(1).LineWidth = 0.1;  h.edge(2).LineWidth = 0.1;
            handles(end+1) = h.mainLine; hold all; %#ok<SAGROW>
        end
    end
%     set(gca, 'XTick', find(ismember(params.xtime,-3e6:0.5e6:3e6)), 'XTickLabels', params.xtime(ismember(params.xtime,(-3e6:0.5e6:3e6)))/1e6,'FontSize', 20)
    set(gca, 'XTick', -3e6:0.5e6:3e6, 'XTickLabels',(-3e6:0.5e6:3e6)/1e6,'FontSize', 20)
    set(gca, 'YTick', [0 0.5 1 1.5], 'FontSize', 20)
    xlim([-0.5e6 1.5e6]);
    ylim([-0.2 1.7])
    
    ylabel('Z-scored firing rate','FontSize', 20)
    xlabel('Time (s)','FontSize', 20)
    legend(handles,params.labels_ztrials); legend boxoff
    
    figure; set(gcf,'units','normalized','Position',[0.65 0.3 0.15 0.23],'color','w'); hold all;
    datatoplot = [];
    datatoplot(1,:) = nanmean(nanmean(zmat(idx_all,timeidx,[1 2]),2),3);
    datatoplot(2,:) = nanmean(nanmean(zmat(idx_all,timeidx,[3 4]),2),3);
    
    errorbar(1,nanmean(datatoplot(1,:),2),nanstd(datatoplot(1,:),[],2) / sqrt(sum(idx_all)),'.','LineWidth',5,'Color',params.colors_ztrials{1});
    errorbar(2,nanmean(datatoplot(2,:),2),nanstd(datatoplot(2,:),[],2) / sqrt(sum(idx_all)),'.','LineWidth',5,'Color',params.colors_ztrials{3});
    
    Y           = [datatoplot(1,:)'; datatoplot(2,:)'];
    X_sal       = [ones(sum(idx_exp),1); ones(sum(idx_exp),1)*2];

    tbl             = table(Y,X_sal,[G_mou(idx_exp); G_mou(idx_exp)],'VariableNames',{'Activity','Saliency','Mouse'}); %Create table for mixed model

    lme             = fitlme(tbl,'Activity~Saliency+(1|Mouse)'); %construct linear mixed effects model with fixed effect of cohort
    stats           = dataset2table(anova(lme,'DFMethod','Satterthwaite')); %Perform ANOVA on model and output as matrix
    fprintf('%s: ',params.ExperimentLabels{iExp})
    fprintf('F(%d,%2.0f)=%1.2f, p=%1.2e; ',stats{2,3},stats{2,4},stats{2,2},stats{2,5})
    sigstar([1 2],stats{2,5})
    ylim([0 0.8])
end
fprintf('Linear Mixed Model ANOVA\n')

%% Make figure of the mean split by PV and PYR:
figure; set(gcf,'units','normalized','Position',[0.2 0.3 0.72 0.38],'color','w')

params.colors_types = {[0.8 0 0] [0 0 0.8]};

params.lines_splits = {'--' '-'};

% suptitle('PV Pyr split')
for iType = 1:2
    iExp = [2 3];
    idx_exp = ismember(spikeData.session_ID,sessionData.session_ID(ismember(sessionData.Experiment,params.Experiments(iExp)))) & spikeData.celltype==iType;
    subplot(1,2,iType)
    sum(idx_exp)
    
    handles = [];
    for iSplit = 3:4 %1:2%length(splits)
        meantoplot = nanmean(zmat(idx_exp,:,iSplit),1);
        errortoplot = nanstd(zmat(idx_exp,:,iSplit),1)/sqrt(size(zmat(idx_exp,:,iSplit),1));
        h = shadedErrorBar(params.xtime,meantoplot,errortoplot,{params.lines_splits{iSplit-2},'markerfacecolor',params.colors_types{iType},'LineWidth',3},0);
        h.mainLine.Color = params.colors_types{iType};    h.patch.FaceColor = params.colors_types{iType};
        delete(h.edge(1)); delete(h.edge(2));
        handles(end+1) = h.mainLine; hold all; %#ok<SAGROW>
    end
   set(gca, 'XTick', -3e6:0.5e6:3e6, 'XTickLabels',(-3e6:0.5e6:3e6)/1e6,'FontSize', 20)
    set(gca, 'YTick', [0 0.5 1 1.5], 'FontSize', 20)
    xlim([-0.5e6 1.5e6]);
    ylim([-0.2 1.9])
    
    ylabel('Z-scored firing rate','FontSize', 20)
    xlabel('Time (s)','FontSize', 20)
    legend(handles,{'Hit' 'Miss'}); legend boxoff
end
