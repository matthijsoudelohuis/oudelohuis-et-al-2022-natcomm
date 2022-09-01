%% This script analyzes the results from the neuro GLM
% This approach fits a generalized linear model to the firing rate based on 
% sensory, decision, task and state variables. The fitting has happened before 
% at a computer cluster. This scripts loads the result of the fitting procedure,
% displays examples and computes measures of data fit quality and variance explained
% Oude Lohuis et al. 2022 Nat Comms

%%
startover

%% Load dataset:
load('E:\Matlab\MOL_Analysis\AOudeLohuisetal_2022_NatComms\3GLM\Data3_1.mat')

%% Output figure parameters
params.nMaxTrials           = max(trialData.trialNum);
params.nPredictors          = size(output.x,3);

params.Experiments          = {'ChangeDetectionConflictDecor' 'VisOnlyTwolevels' 'ChangeDetectionConflict' }; %Which versions of the task to load data from
params.ExperimentLabels     = {'NE' 'UST' 'MST'}; %Which versions of the task to load data from

params                      = MOL_getColors_CHDET(params);

params.posthoctest          = 'bonferroni';
params.alpha                = 0.01;

params.savedir              = 'E:\Documents\PhD\Figures\Project CHDET\Results - 2nd Bump\3 neuroGLM';

%% Compute variance explained over all individual trials:
idx         = params.xtime>0 & params.xtime<2e6;
idx         = repmat(idx,1,params.nMaxTrials);
for iNeuron = 1:params.nNeurons %loop over neurons
    for iM = 1:params.nModels %loop over models
        %Compute variance explained by the model for this neuron (across certain time range)
        output.var_expl_v1(iNeuron,iM)     = 1 - nanvar(output.y(iNeuron,idx)' - squeeze(output.y_hat(iNeuron,iM,idx))) / nanvar(output.y(iNeuron,idx));
    end
end

output.var_expl_v1(output.var_expl_v1<0) = 0;

%% Compute variance explained over trialtype-averaged activity (Runyan et al. Harvey lab custom)
idx         = params.xtime>0 & params.xtime<2e6;
splitidxs= {};
for iNeuron = 1:params.nNeurons %loop over neurons
    [temptrialData]     = MOL_getTempPerSes(spikeData.session_ID(iNeuron),trialData);
    nTotalSpikeBins     = sum(~isnan(output.y(iNeuron,:)));
    nTrials             = nTotalSpikeBins/params.nTimebins;
    hist_mat            = reshape(output.y(iNeuron,1:nTotalSpikeBins),params.nTimebins,nTrials);
    
    splitidxs{1} = strcmp(temptrialData.trialType,'X') & temptrialData.visualOriChangeNorm==3 & temptrialData.vecResponse==3;
    splitidxs{2} = strcmp(temptrialData.trialType,'X') & temptrialData.visualOriChangeNorm==3 & temptrialData.vecResponse==2;
    splitidxs{3} = strcmp(temptrialData.trialType,'Y') & temptrialData.audioFreqChangeNorm==3 & temptrialData.vecResponse==3;
    splitidxs{4} = strcmp(temptrialData.trialType,'Y') & temptrialData.audioFreqChangeNorm==3 & temptrialData.vecResponse==1;
    splitidxs{5} = strcmp(temptrialData.trialType,'P') & temptrialData.vecResponse==3;
    splitidxs{6} = strcmp(temptrialData.trialType,'P') & ismember(temptrialData.vecResponse, [1 2]);

    hist_mat_mean = [];
    for iSplit = 1:length(splitidxs)
        hist_mat_mean = [hist_mat_mean; nanmean(hist_mat(idx,splitidxs{iSplit}),2)]; %#ok<*AGROW>
    end
    
    for iM = 1:params.nModels
        hist_model          = reshape(output.y_hat(iNeuron,iM,1:nTotalSpikeBins),params.nTimebins,nTrials);
        
        hist_model_mean = [];
        for iSplit = 1:length(splitidxs)
            hist_model_mean = [hist_model_mean; nanmean(hist_model(idx,splitidxs{iSplit}),2)];
        end
        output.var_expl_v2(iNeuron,iM)     = 1 - nanvar(hist_mat_mean - hist_model_mean) / nanvar(hist_mat_mean);
    end
end

output.var_expl_v2(output.var_expl_v2<0) = 0; %set var explained below 0 to 0

%% Create figure PSTH of some example cells that have a high fit quality:
% %take some examples from each experiment type:
% temp                        = output.var_expl_v2(:,2) - output.var_expl_v2(:,1);
% params.example_cell_IDs     = []; %init var to store some cell ids of example neurons
% for iExp = 1:length(params.Experiments)
%     expidx = ismember(spikeData.session_ID,sessionData.session_ID(strcmp(sessionData.Experiment,params.Experiments(iExp))));
%     params.example_cell_IDs = [params.example_cell_IDs; spikeData.cell_ID(temp > prctile(temp(expidx),95,1) & expidx)]; %with high expl variance
% end

%Set the examples fixed:
params.example_cell_IDs = {'10092019030831139' '20282019121311450' '20112018081011056'};

for iNeuron = find(ismember(spikeData.cell_ID,params.example_cell_IDs))' %Show for example neuron i:
    
    [temptrialData]     = MOL_getTempPerSes(spikeData.session_ID(iNeuron),trialData);
    
    nTotalSpikeBins     = sum(~isnan(output.y(iNeuron,:)));
    nTrials             = nTotalSpikeBins/params.nTimebins;
    hist_mat            = reshape(output.y(iNeuron,1:nTotalSpikeBins),params.nTimebins,nTrials);
    
    figure; hold all;
    set(gcf,'units','normalized','Position',[0.2 0.1 0.3 0.4],'color','w')
    title(sprintf('%s, EV1: %2.0f EV2: %2.0f',spikeData.cell_ID{iNeuron},output.var_expl_v1(iNeuron,2)*100,output.var_expl_v2(iNeuron,2)*100),'FontSize',10)

    colors      = {[0.4 0.4 1]      [0 0 1]         [1 0.4 0.4]     [1 0 0]     [0.7 0.7 0.7]       [0.3 0.3 0.3]};
    labels      = {'Visual Miss'    'Visual Hit'    'Audio Miss'    'Audio Hit' 'Correct Rejection' 'False Alarm'};
    %Show for these trial types the avg firing rate and model estimated firing rate:
    splitidxs= {};
    splitidxs{1} = strcmp(temptrialData.trialType,'X') & temptrialData.visualOriChangeNorm==3 & temptrialData.vecResponse==3;
    splitidxs{2} = strcmp(temptrialData.trialType,'X') & temptrialData.visualOriChangeNorm==3 & temptrialData.vecResponse==2;
    splitidxs{3} = strcmp(temptrialData.trialType,'Y') & temptrialData.audioFreqChangeNorm==3 & temptrialData.vecResponse==3;
    splitidxs{4} = strcmp(temptrialData.trialType,'Y') & temptrialData.audioFreqChangeNorm==3 & temptrialData.vecResponse==1;
    splitidxs{5} = strcmp(temptrialData.trialType,'P') & temptrialData.vecResponse==3;
    splitidxs{6} = strcmp(temptrialData.trialType,'P') & ismember(temptrialData.vecResponse, [1 2]);

    iM = 2; %show only for full model:
    hist_model          = reshape(output.y_hat(iNeuron,iM,1:nTotalSpikeBins),params.nTimebins,nTrials);
    CurrLineHandles = [];
    for iSplit = 1:length(splitidxs)
        
        %Get the mean over trials::
        hist_mat_mean           = mean(hist_mat(:,splitidxs{iSplit}),2);
        hist_mat_model_mean     = mean(hist_model(:,splitidxs{iSplit}),2);
        
        CurrLineHandles(end+1) = plot(params.xtime,hist_mat_mean,'-','LineWidth',2,'Color',colors{iSplit}); %#ok<SAGROW>
        plot(params.xtime,hist_mat_model_mean,':','LineWidth',2,'Color',colors{iSplit});
        
    end
    
    %Figure makeup:
    xlabel('Time (sec)','FontSize',15)    %Label x-axis
    ylabel('Hz','FontSize',15)            %Label y-axis
    set(gca, 'XTick', params.xtime(1:500e3/params.binsize:params.nTimebins), 'XTickLabels', params.xtime(1:500e3/params.binsize:params.nTimebins)/1e6,'FontSize', 20)
    set(gca, 'FontSize', 20);
    xlim([params.xtime(1) params.xtime(end)]);
    legend(CurrLineHandles,labels,'FontSize',15);
    legend(gca,'boxoff');
    
    export_fig(fullfile(params.savedir,sprintf('ExNeuron_ModelFull_ExplVar_%s.eps',spikeData.cell_ID{iNeuron})))
end

%% For statistics:
X_coh           = cell(params.nNeurons,1);
X_coh(ismember(spikeData.session_ID,sessionData.session_ID(ismember(sessionData.Experiment,params.Experiments{1}))))         = params.ExperimentLabels(1);
X_coh(ismember(spikeData.session_ID,sessionData.session_ID(ismember(sessionData.Experiment,params.Experiments{2}))))         = params.ExperimentLabels(2);
X_coh(ismember(spikeData.session_ID,sessionData.session_ID(ismember(sessionData.Experiment,params.Experiments{3}))))         = params.ExperimentLabels(3);

G_mou           = NaN(params.nNeurons,1);
uMice           = unique(sessionData.mousename);
for iMouse = 1:length(uMice)
    G_mou(ismember(spikeData.session_ID,sessionData.session_ID(strcmp(sessionData.mousename,uMice{iMouse})))) = iMouse;
end

%% Figure variance explained:
xpos = [2 4 6];
datatoplot = NaN(6,1000);
for iExp = 1:length(params.Experiments)
    expidx = ismember(spikeData.session_ID,sessionData.session_ID(strcmp(sessionData.Experiment,params.Experiments(iExp))));
    datatoplot((iExp-1)*2+1,1:sum(expidx)) = output.var_expl_v1(expidx,1);
    datatoplot((iExp-1)*2+2,1:sum(expidx)) = output.var_expl_v1(expidx,2);
end
%     violinplot(datatoplot',[],'Width',0.45,'ViolinColor',params.colors_experiments,'ViolinAlpha',1,'ShowData',false);
figure; hold all; set(gcf,'units','normalized','Position',[0.35 0.45 0.15 0.28],'color','w');
% violinplot(datatoplot',[],'Width',0.45,'ViolinColor',params.colors_experiments{iExp},'ViolinAlpha',1,'ShowData',false);
% boxplot(datatoplot',[],'Width',0.45,'Color',params.colors_experiments{iExp});
boxplot(datatoplot','plotstyle','compact','Color',params.colors_experiments{iExp},'boxstyle','filled', 'medianstyle','target','whisker',5);
set(gca,'XTickLabels',params.modelString,'YTick',[0 0.4],'XTick',[2 4 6],'XTickLabels',params.ExperimentLabels)
ylim([0 0.4]);
ylabel('Explained variance')

%Statistical testing between cohorts:
Y               = output.var_expl_v1(:,2);
tbl             = table(Y,X_coh,'VariableNames',{'EV','Cohort'}); %Create table for mixed model
lme             = fitlme(tbl,'EV~Cohort'); %construct linear mixed effects model with fixed effect of temporal window and random intercept for different mice
stats           = dataset2table(anova(lme,'DFMethod','Satterthwaite')); %Perform ANOVA on model and output as matrix
fprintf('Fixed effect of cohort on EV: (Linear Mixed Model ANOVA)\n')
fprintf('(F(%d,%2.0f) = %1.2f, p=%1.2f, ANOVA)\n',stats{2,3},stats{2,4},stats{2,2},stats{2,5})

contrasts           = {[0 1 0] [0 1 -1] [0 0 1]};
contrastlabels      = {'MST vs NE' 'NE vs UST' 'MST vs UST'};

fprintf('Posthoc comparison:\n')
for iC = 1:3
    [p,F,DF1,DF2] = coefTest(lme,contrasts{iC});
    fprintf('(%s):F(%d,%d)=%1.2f, p=%1.3f\n',contrastlabels{iC},DF1,DF2,F,p)
end

tempfile = fullfile('SourceData_FigS4a_GLM_EV_v1_Cohorts.xlsx');
writetable(tbl,tempfile)

%% Figure variance explained:
datatoplot = NaN(6,1000);
for iExp = 1:length(params.Experiments)
    expidx = ismember(spikeData.session_ID,sessionData.session_ID(strcmp(sessionData.Experiment,params.Experiments(iExp))));
    datatoplot((iExp-1)*2+1,1:sum(expidx)) = output.var_expl_v2(expidx,1);
    datatoplot((iExp-1)*2+2,1:sum(expidx)) = output.var_expl_v2(expidx,2);
end
%     violinplot(datatoplot',[],'Width',0.45,'ViolinColor',params.colors_experiments,'ViolinAlpha',1,'ShowData',false);
figure; hold all; set(gcf,'units','normalized','Position',[0.65 0.45 0.15 0.28],'color','w');

boxplot(datatoplot','plotstyle','compact','Color',params.colors_experiments{iExp},'boxstyle','filled', 'medianstyle','target','whisker',5);
set(gca,'XTickLabels',params.modelString,'YTick',[0 0.5 1],'XTick',[2 4 6],'XTickLabels',params.ExperimentLabels)
ylim([0 1])
ylabel('Explained variance')

%Statistical testing between cohorts:
Y               = output.var_expl_v2(:,2);
tbl             = table(Y,X_coh,'VariableNames',{'EV','Cohort'}); %Create table for mixed model
lme             = fitlme(tbl,'EV~Cohort'); %construct linear mixed effects model with fixed effect of temporal window and random intercept for different mice
stats           = dataset2table(anova(lme,'DFMethod','Satterthwaite')); %Perform ANOVA on model and output as matrix
fprintf('Fixed effect of cohort on EV: (Linear Mixed Model ANOVA)\n')
fprintf('(F(%d,%2.0f) = %1.2f, p=%1.3f, ANOVA)\n',stats{2,3},stats{2,4},stats{2,2},stats{2,5})

contrasts           = {[0 1 0] [0 1 -1] [0 0 1]};
contrastlabels      = {'MST vs NE' 'NE vs UST' 'MST vs UST'};

fprintf('Posthoc comparison:\n')
for iC = 1:3
    [p,F,DF1,DF2] = coefTest(lme,contrasts{iC});
    fprintf('(%s):F(%d,%d)=%1.2f, p=%1.3f\n',contrastlabels{iC},DF1,DF2,F,p)
end

tempfile = fullfile('SourceData_FigS4b_GLM_EV_v2_Cohorts.xlsx');
writetable(tbl,tempfile)

%% Figure: var explained as a function of time for subselection of variables:
var_expl_overtime       = NaN(params.nNeurons,params.nTimebins,params.nSplits);
for iS = 1:params.nSplits
    for iNeuron = 1:params.nNeurons %loop over neurons
        nTotalSpikeBins                 = sum(~isnan(output.y(iNeuron,:)));
        nTrials                         = nTotalSpikeBins/params.nTimebins;
        hist_mat                        = reshape(output.y(iNeuron,1:nTotalSpikeBins),params.nTimebins,nTrials);
        hist_mat_hat                    = reshape(output.y_hat_split(iNeuron,iS,1:nTotalSpikeBins),params.nTimebins,nTrials);
        var_expl_overtime(iNeuron,:,iS) = 1 - nanvar(hist_mat - hist_mat_hat,[],2) ./ nanvar(hist_mat,[],2);
    end
end

%Set to NaN impossible values (more than 100% explained, or 80% more variance due to model:
var_expl_overtime(var_expl_overtime>1 | var_expl_overtime<-0.8) = NaN;

%% Make the figure:
idx_vars = [2 5 4 6];

figure; set(gcf,'units','normalized','Position',[0.05 0.5 0.9 0.3],'color','w'); hold all;
for iS = 1:length(idx_vars)
    iVar = idx_vars(iS);
    subplot(1,length(idx_vars),iS); hold all;
    handles = [];
    
    datatoplot = NaN(1000,params.nTimebins,params.nExperiments); 
    for iExp = 1:length(params.Experiments)
        expidx = ismember(spikeData.session_ID,sessionData.session_ID(strcmp(sessionData.Experiment,params.Experiments(iExp))));
        if iVar==5 && iExp==2; expidx = ismember(spikeData.session_ID,sessionData.session_ID(strcmp(sessionData.Experiment,params.Experiments(3)))); 
        elseif iVar==5 && iExp==3; expidx = ismember(spikeData.session_ID,sessionData.session_ID(strcmp(sessionData.Experiment,params.Experiments(2))));
        end
        
        datatoplot(1:sum(expidx),:,iExp)    = var_expl_overtime(expidx,:,iVar); %get for each cohort the variance explained over time in the datatoplot variable
        meantoplot                          = squeeze(nanmean(datatoplot(:,:,iExp),1)); %take mean and trim dimensions:
        errortoplot                         = squeeze(nanstd(datatoplot(:,:,iExp),1)) / sqrt(sum(expidx));
        h = shadedErrorBar(params.xtime*1e-6,meantoplot,errortoplot,{'-k','markerfacecolor',params.colors_experiments{iExp},'LineWidth',2},0);
        %line make up:
        h.mainLine.Color = params.colors_experiments{iExp};    h.patch.FaceColor = params.colors_experiments{iExp};
        delete(h.edge(1)); delete(h.edge(2));
        handles(end+1) = h.mainLine;
    end
    title(sprintf('%s',strrep(params.varsplits{iVar},'var_',' '))) %set title string of variable

    %Significance is not in the main figure:
    pthr = 0.05 / params.nTimebins;
    for iBin = 1:params.nTimebins
        [~,h_NE_UST(iBin)]      = ranksum(datatoplot(:,iBin,1),datatoplot(:,iBin,2),'alpha',pthr); %#ok<*SAGROW>
        [~,h_NE_MST(iBin)]      = ranksum(datatoplot(:,iBin,1),datatoplot(:,iBin,3),'alpha',pthr);
        [~,h_UST_MST(iBin)]     = ranksum(datatoplot(:,iBin,2),datatoplot(:,iBin,3),'alpha',pthr);
    end
    plot(params.xtime(h_NE_UST)*1e-6,repmat(0.048,sum(h_NE_UST),1),'ks','LineWidth',1.5,'MarkerFaceColor','k')
    plot(params.xtime(h_NE_MST)*1e-6,repmat(0.046,sum(h_NE_MST),1),'ks','LineWidth',1.5,'MarkerFaceColor','k')
    plot(params.xtime(h_UST_MST)*1e-6,repmat(0.044,sum(h_UST_MST),1),'ks','LineWidth',1.5,'MarkerFaceColor','k')
    
    xlim([-0.5 1.5])
    ylim([0 0.05])
    set(gca,'XTick',[-0.5 0 0.5 1 1.5],'YTick',[0 0.05])
    plot([0e6 0e6],[0 1],'--','Color',[0.5 0.5 0.5],'LineWidth',2)
    if iS==1
        legend(handles,params.ExperimentLabels); legend boxoff
    end
end

%% Repeat onset latency analysis as with AUC:

idx_vars            = [2 5]; %focus on visual and hit/miss coding

idx_resp            = params.xtime>0e6 & params.xtime<1e6; %only during response window
idx_bsl             = params.xtime<0e6; %use baseline for threshold computation

%Define threshold as baseline + 2 standard deviations
tempbsl             = var_expl_overtime(:,idx_bsl,idx_vars);
params.thr          = nanmean(tempbsl(:))+2*nanstd(tempbsl(:));

%init output:
var_expl_thr        = NaN(params.nNeurons,params.nSplits);

for iNeuron = 1:params.nNeurons %loop over neurons
    for iS = 1:params.nSplits %For each variable define first moment variance explained exceeds threshold
        temp = params.xtime(find(var_expl_overtime(iNeuron,:,iS)>params.thr & idx_resp,1));
        if temp %if exceeds then store
            var_expl_thr(iNeuron,iS)        =  temp;
        end
    end
end

figure; set(gcf,'units','normalized','Position',[0.1 0.4 0.35 0.45],'color','w'); hold all;
exporder = [1 3 2];
for iExp = 1:3 %show for each cohort the center and distr of earliest onset lat of encoding:
    idx_exp             = ismember(spikeData.session_ID,sessionData.session_ID(strcmp(sessionData.Experiment,params.Experiments(exporder(iExp)))));
    
    meantoplot          = squeeze(nanmedian(var_expl_thr(idx_exp,:),1));
    errortoplot         = squeeze(nanstd(var_expl_thr(idx_exp,:),[],1)) / sqrt(sum(idx_exp));

    errorbar((1:length(idx_vars))-0.3+iExp*0.2,meantoplot(idx_vars),errortoplot(idx_vars),'.','color',params.colors_experiments{iExp},'LineWidth',3,'MarkerSize',50);
    view([90 -90])
    grid on;
end

% Figure make up:
set(gca,'XTick',1:length(idx_vars),'XTickLabels',params.varsplits(idx_vars),'FontSize', 25)
set(gca,'Position',[0.25 0.2 0.7 0.7])
ylim([0 0.45e6])
set(gca,'YTick',0:0.1e6:0.4e6,'YTickLabels',(0:0.1e6:0.4e6)*1e-3,'Fontsize',15);
ylabel('Earliest time of significant fraction neurons (ms)','Fontsize',15);
legend(params.ExperimentLabels([1 2 3]),'Location','northeast'); legend boxoff

%Multi-level statistics: 
X_coh           = NaN(params.nNeurons,1);
X_coh(ismember(spikeData.session_ID,sessionData.session_ID(ismember(sessionData.Experiment,params.Experiments{1}))))         = 1;
X_coh(ismember(spikeData.session_ID,sessionData.session_ID(ismember(sessionData.Experiment,params.Experiments{2}))))         = 2;
X_coh(ismember(spikeData.session_ID,sessionData.session_ID(ismember(sessionData.Experiment,params.Experiments{3}))))         = 3;

Y_lat           = var_expl_thr(:,5);
idx             = ismember(X_coh,[2 3]); %only compare UST and MST
tbl             = table(Y_lat(idx),X_coh(idx),G_mou(idx),'VariableNames',{'EV','Cohort','Mouse'}); %Create table for mixed model

lme             = fitlme(tbl,'EV~Cohort'); %construct linear mixed effects model with fixed effect of temporal window and random intercept for different mice
stats           = dataset2table(anova(lme,'DFMethod','Satterthwaite')); %Perform ANOVA on model and output as matrix
fprintf('Onset UST vs MST: (Linear Mixed Model)\n')
fprintf('(F(%d,%2.0f) = %1.2f, p=%1.2f, ANOVA)\n',stats{2,3},stats{2,4},stats{2,2},stats{2,5})

tempfile = fullfile('SourceData_FigS4e_GLM_EVonset_Cohorts.xlsx');
writetable(tbl,tempfile)

%% Compute variance explained over all individual trials for subselection of variables:
idx         = params.xtime>0e6 & params.xtime<=0.2e6;
idx         = repmat(idx,1,params.nMaxTrials);

var_expl_splits_early       = NaN(params.nNeurons,params.nSplits);
for iS = 1:params.nSplits
    for iNeuron = 1:params.nNeurons %loop over neurons
        %Compute variance explained by this split (subselection of predictors) for this neuron (across certain time range)
        var_expl_splits_early(iNeuron,iS)     = 1 - nanvar(output.y(iNeuron,idx)' - squeeze(output.y_hat_split(iNeuron,iS,idx))) / nanvar(output.y(iNeuron,idx));
    end
end

%Exclude infinite value occuring in a neuron;
var_expl_splits_early(var_expl_splits_early>1 | var_expl_splits_early<-1) = NaN;

%% Make figure;
idx_vars        = [2 4 5 6];
nVars           = length(idx_vars);

%Init plotvar:
datatoplot          = NaN(500,params.nExperiments,params.nSplits);
for iExp = 1:params.nExperiments %For each experiment copute the mean and error by averaging expl var over neurons
    for iS = 1:params.nSplits
        expidx                              = ismember(spikeData.session_ID,sessionData.session_ID(strcmp(sessionData.Experiment,params.Experiments(iExp))));
        if iS==5 && iExp==2
            expidx = ismember(spikeData.session_ID,sessionData.session_ID(strcmp(sessionData.Experiment,params.Experiments(3))));
        elseif iS==5 && iExp==3
            expidx = ismember(spikeData.session_ID,sessionData.session_ID(strcmp(sessionData.Experiment,params.Experiments(2))));
        end
        datatoplot(1:sum(expidx),iExp,iS)           = var_expl_splits_early(expidx,iS);
    end
end

clear h;
figure; hold all; set(gcf,'units','normalized','Position',[0.05 0.5 0.2 0.2],'color','w');
title('Variance explained - 0-200ms')

for i = 1:nVars
    xpos = [0.8 1 1.2] + (i-1);
    
    h = boxplot(datatoplot(:,:,idx_vars(i)),'positions',xpos,'outliersize',0.0001,'widths',0.16);
    for j=1:3
        set(h(1,j),'Color',params.colors_experiments{j},'LineStyle','-','LineWidth',1);
        set(h(2,j),'Color',params.colors_experiments{j},'LineStyle','-','LineWidth',1);
        set(h(3,j),'Color',params.colors_experiments{j});
        set(h(4,j),'Color',params.colors_experiments{j});
        set(h(5,j),'Color',params.colors_experiments{j},'LineWidth',1);
        set(h(6,j),'Color',params.colors_experiments{j},'LineWidth',1);
        delete(h(7,j));
    end
    
    %Statistical testing:
    Y               = var_expl_splits_early(:,idx_vars(i));
    tbl             = table(Y,X_coh,'VariableNames',{'EV','Cohort'}); %Create table for mixed model

    lme             = fitlme(tbl,'EV~Cohort'); %construct linear mixed effects model with fixed effect of cohort
    stats           = dataset2table(anova(lme,'DFMethod','Satterthwaite')); %Perform ANOVA on model and output as matrix
    fprintf('Var expl early per cohort (%s): (Linear Mixed Model)\n',params.varsplits{idx_vars(i)})
    fprintf('(F(%d,%2.0f) = %1.2f, p=%1.2f, ANOVA)\n',stats{2,3},stats{2,4},stats{2,2},stats{2,5})

    ordering = [3 1 2]; %annoying reordering of the groups by lme (MST is taken as first reference category)
    if stats{2,5}<0.05
        contrasts           = {[0 1 0] [0 1 -1] [0 0 1]};
        contrastlabels      = {'MST vs NE' 'NE vs UST' 'MST vs UST'};
        
        fprintf('Posthoc comparison:\n')
        for iC = 1:3
            [p,F,DF1,DF2] = coefTest(lme,contrasts{iC});
            fprintf('(%s):F(%d,%d)=%1.2f, p=%1.3f\n',contrastlabels{iC},DF1,DF2,F,p)
            if p<0.05
               if sum(contrasts{iC})==1
                   contrasts{iC}(1) = 1; %put MST as position for sigstar in contrast var (not used for other purposes)
               end
                sigstar(xpos(ordering(contrasts{iC}~=0)),p) %use sigstar function to identify
            end
        end
    end
    
    tempfile = fullfile(sprintf('SourceData_FigS4c_GLM_EV_Early_%s.xlsx',params.varsplits{idx_vars(i)}));
    writetable(tbl,tempfile)
    
end

set(gca,'XTick',1:nVars,'XTickLabel',strrep(params.varsplits(idx_vars),'_',' '),'XTickLabelRotation',0,'Fontsize',10,'YTick',[0 0.025 0.05 0.075 0.1]);
% legend(h,params.ExperimentLabels); legend boxoff;
ylabel('Explained variance')
% ylim([0 0.06])
ylim([-0.01 0.1])
xlim([0.5 nVars+0.5])

% % saveas(gcf,fullfile(params.savedir,'V1_neuroGLM_VarExpl_perVar_Bar_Early.eps'))
export_fig(fullfile(params.savedir,'V1_neuroGLM_VarExpl_perVar_Box_Early.eps'),gcf)

%% Compute variance explained over all individual trials for subselection of variables:
idx         = params.xtime>0.2e6 & params.xtime<=1.5e6;
idx         = repmat(idx,1,params.nMaxTrials);

var_expl_splits_late       = NaN(params.nNeurons,params.nSplits);
for iS = 1:params.nSplits
    for iNeuron = 1:params.nNeurons %loop over neurons
        %Compute variance explained by this split (subselection of predictors) for this neuron (across certain time range)
        var_expl_splits_late(iNeuron,iS)     = 1 - nanvar(output.y(iNeuron,idx)' - squeeze(output.y_hat_split(iNeuron,iS,idx))) / nanvar(output.y(iNeuron,idx));
    end
end

%Exclude infinite value occuring in a neuron;
var_expl_splits_late(var_expl_splits_late>1 | var_expl_splits_late<-1) = NaN;

%% Make figure;
idx_vars        = [2 4 5 6];
nVars           = length(idx_vars);


%Init plotvar:
datatoplot          = NaN(500,params.nExperiments,params.nSplits);
for iExp = 1:params.nExperiments %For each experiment copute the mean and error by averaging expl var over neurons
    for iS = 1:params.nSplits
        expidx                              = ismember(spikeData.session_ID,sessionData.session_ID(strcmp(sessionData.Experiment,params.Experiments(iExp))));
        if iS==5 && iExp==2
            expidx = ismember(spikeData.session_ID,sessionData.session_ID(strcmp(sessionData.Experiment,params.Experiments(3))));
        elseif iS==5 && iExp==3
            expidx = ismember(spikeData.session_ID,sessionData.session_ID(strcmp(sessionData.Experiment,params.Experiments(2))));
        end
        datatoplot(1:sum(expidx),iExp,iS)           = var_expl_splits_late(expidx,iS);
    end
end


clear h;
figure; hold all; set(gcf,'units','normalized','Position',[0.05 0.5 0.2 0.2],'color','w');
title('Variance explained - 200-1000ms')
for i = 1:nVars
    xpos = [0.8 1 1.2] + (i-1);
%     xpos = [0.7 0.9 1.1 1.3] + (i-1);


    h = boxplot(datatoplot(:,:,idx_vars(i)),'positions',xpos,'outliersize',0.0001,'widths',0.16);
    for j=1:3
        set(h(1,j),'Color',params.colors_experiments{j},'LineStyle','-','LineWidth',1);
        set(h(2,j),'Color',params.colors_experiments{j},'LineStyle','-','LineWidth',1);
        set(h(3,j),'Color',params.colors_experiments{j});
        set(h(4,j),'Color',params.colors_experiments{j});
        set(h(5,j),'Color',params.colors_experiments{j},'LineWidth',1);
        set(h(6,j),'Color',params.colors_experiments{j},'LineWidth',1);
        delete(h(7,j));
    end

    %Statistical testing:
    Y               = var_expl_splits_late(:,idx_vars(i));
    tbl             = table(Y,X_coh,'VariableNames',{'EV','Cohort'}); %Create table for mixed model

    lme             = fitlme(tbl,'EV~Cohort'); %construct linear mixed effects model with fixed effect of cohort
    stats           = dataset2table(anova(lme,'DFMethod','Satterthwaite')); %Perform ANOVA on model and output as matrix
    fprintf('Var expl early per cohort (%s): (Linear Mixed Model)\n',params.varsplits{idx_vars(i)})
    fprintf('(F(%d,%2.0f) = %1.2f, p=%1.2f, ANOVA)\n',stats{2,3},stats{2,4},stats{2,2},stats{2,5})

    ordering = [3 1 2]; %annoying reordering of the groups by lme (MST is taken as first reference category)
    if stats{2,5}<0.05
        contrasts           = {[0 1 0] [0 1 -1] [0 0 1]};
        contrastlabels      = {'MST vs NE' 'NE vs UST' 'MST vs UST'};
        
        fprintf('Posthoc comparison:\n')
        for iC = 1:3
            [p,F,DF1,DF2] = coefTest(lme,contrasts{iC});
            fprintf('(%s):F(%d,%d)=%1.2f, p=%1.3f\n',contrastlabels{iC},DF1,DF2,F,p)
            if p<0.05
               if sum(contrasts{iC})==1
                   contrasts{iC}(1) = 1; %put MST as position for sigstar in contrast var (not used for other purposes)
               end
                sigstar(xpos(ordering(contrasts{iC}~=0)),p) %use sigstar function to identify
            end
        end
    end
    
    tempfile = fullfile(sprintf('SourceData_FigS4d_GLM_EV_Late_%s.xlsx',params.varsplits{idx_vars(i)}));
    writetable(tbl,tempfile)
end

set(gca,'XTick',1:nVars,'XTickLabel',strrep(params.varsplits(idx_vars),'_',' '),'XTickLabelRotation',0,'Fontsize',10,'YTick',[0 0.025 0.05 0.075 0.1]);
% legend(h,params.ExperimentLabels); legend boxoff;
ylabel('Explained variance')
ylim([-0.01 0.1])
xlim([0.5 nVars+0.5])

export_fig(fullfile(params.savedir,'V1_neuroGLM_VarExpl_perVar_Box_Late.eps'),gcf)

