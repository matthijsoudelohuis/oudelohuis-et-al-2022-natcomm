%% This script analyzes lick related activity
% This script aligns V1 neural activity to licking onset to relate spiking 
% activity to the timing of goal-directed movement
% Supplementary Figure to:
% Oude Lohuis et al. 2022 Nat Comms

%% Reset all
startover

%% Parameter settings:
params                  = params_histresponse(); %Parameters for PSTH (All time is in microseconds)

params.Experiments      = {'ChangeDetectionConflictDecor' 'VisOnlyTwolevels' 'ChangeDetectionConflict' }; %Which versions of the task to load data from
params.ExperimentLabels = {'NE' 'UST' 'MST'}; %Labels for the different experiments

params.AlignOn          = 'stimChange';      %On which timestamp to align as t=0

params.minTrialCond     = 3;
params.cscale           = [-0.5 2.2];

params.area             = 'V1';


params                  = MOL_getColors_CHDET(params);
params.savedir          = 'E:\Documents\PhD\Figures\Project CHDET\Results - 2nd Bump\2 Neural\Zmat_AlignLick';

%% Load dataset:
load('E:\Matlab\MOL_Analysis\AOudeLohuisetal_2022_NatComms\2Neural\Data2_1.mat')
fprintf('Dataset: %d sessions, %d trials, %d neurons\n',length(sessionData.session_ID),length(trialData.session_ID),length(spikeData.session_ID));

%% Set parameters:
params                      = params_histresponse_coding(params);
params.t_pre                = -3e6;
params.t_post               = 3e6;
params.twin_baseline_start  = -3e6;
params.twin_baseline_stop   = -2e6;

%Construct bin edges and time axis
params.edges                = [params.t_pre:params.binsize:params.t_post] - params.binsize/2;                    %#ok<NBRAK> Define edges of the bins
params.xtime                = params.t_pre:params.binsize:params.t_post-params.binsize;                        %Define time axis
params.nTimebins            = length(params.xtime); %number of time bins

params.zscore               = 1;

params.AlignOn              = 'firstLick';      %On which timestamp to align as t=0

params.nSplits              = 9;
params.trialcolors          = {[0.7 0 0] [0.2 0.2 1] [0.7 0.5 0.5] [1 0.2 0.2] [0 0 0.7] [0.5 0.5 0.7] [1 0.5 0.5] [0.5 0.5 1] [0 0 0] };
params.triallines           = {'-' '-' ':' '-' '-' ':' '-' '-' ':' };
params.triallabels          = {'Audio Hit' 'Audio Error' 'Audio Miss' 'Visual Error' 'Visual Hit' 'Visual Miss' 'FA Aud' 'FA Vis' 'CorrRej' };

%Make new field with first lick latency to align on. For probe trials make surrogate timestamp based on random latency
nMiss                       = sum(isnan(trialData.responseLatency));
respmean                    = nanmean(trialData.responseLatency);
respstd                     = nanstd(trialData.responseLatency);
%Create field 'firstLick':
trialData.responseLatency(isnan(trialData.responseLatency)) = randn(nMiss,1)*respstd/2+respmean;
trialData.firstLick         = trialData.stimChange + trialData.responseLatency;

%% Main loop to get psth matrix:

nNeurons                = length(spikeData.ts);
lastsesid               = []; %memory var to keep track of whether neuron comes from different session and design matrix needs to be reconstructed
fprintf('Computing average Z-scored response for neuron        \n');
zmat                = NaN(nNeurons,params.nTimebins,params.nSplits);

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
    splits{end+1}       = ismember(temptrialData.trialType,{'Y'}) & temptrialData.vecResponse==1; %#ok<*SAGROW>
    splits{end+1}       = ismember(temptrialData.trialType,{'Y'}) & temptrialData.vecResponse==2;
    splits{end+1}       = ismember(temptrialData.trialType,{'Y'}) & temptrialData.vecResponse==3;
    splits{end+1}       = ismember(temptrialData.trialType,{'X'}) & temptrialData.vecResponse==1;
    splits{end+1}       = ismember(temptrialData.trialType,{'X'}) & temptrialData.vecResponse==2;
    splits{end+1}       = ismember(temptrialData.trialType,{'X'}) & temptrialData.vecResponse==3;
    splits{end+1}       = ismember(temptrialData.trialType,{'P'}) & temptrialData.vecResponse==1;
    splits{end+1}       = ismember(temptrialData.trialType,{'P'}) & temptrialData.vecResponse==2;
    splits{end+1}       = ismember(temptrialData.trialType,{'P'}) & temptrialData.vecResponse==3;
    
    for iSplit = 1:params.nSplits
        if sum(splits{iSplit})>=params.minTrialCond
            zmat(iNeuron,:,iSplit) = nanmean(hist_mat(splits{iSplit},:),1);
        end
    end
end

%% Sort zmat by averaged activity around lick time and sort by visual and auditory hits
params.twin_resp_start  = -0.3e6;
params.twin_resp_stop   = 0.3e6;

idx                     = [1 5]; %Visual and auditory hits

meanresp                = nanmean(nanmean(zmat(:,params.xtime>params.twin_resp_start & params.xtime<=params.twin_resp_stop,idx),2),3);
[~,sortidx]             = sort(meanresp,1,'descend');
zmat                    = zmat(sortidx,:,:);
spikeFields             = fieldnames(spikeData);
for iF                  = 1:length(spikeFields)%also sort the spikeData to know which cell belongs to which session etcetera
    spikeData.(spikeFields{iF}) = spikeData.(spikeFields{iF})(sortidx,:);
end

%% Make figure:

zmat(isnan(zmat))       = 0;
cmat                    = parula(10000);
temp                    = linspace(params.cscale(1),params.cscale(2),10000);
temp(find(temp>0,1)-1)  = 0;
cmat(temp==0,:)         = [0.9 0.9 0.9]; %make colormap that blanks out trials with no activity

nExperiments = length(params.Experiments);
for iExp = 1:nExperiments
    idx_exp = ismember(spikeData.session_ID,sessionData.session_ID(strcmp(sessionData.Experiment,params.Experiments(iExp))));
    
    figure; set(gcf,'units','normalized','Position',[0.1 0.02 0.65 0.9],'color','w')
%     suptitle(params.ExperimentLabels{iExp})
    for iSplit = 1:params.nSplits
        ax = subplot(3,3,iSplit);
%         pcolor(zmat(idx_exp,:,iSplit)); hold on;
%         shading flat
%         set(gca, 'ydir', 'reverse');
%         set(ax, 'clim', params.cscale);

        temp = zmat(idx_exp,:,iSplit); temp = temp(~all(temp==0,2),:);

        imagesc(temp,params.cscale); hold on;
%         imagesc(zmat(idx_exp,:,iSplit),params.cscale); hold on;
        colormap(cmat)
        plot([find(params.xtime == 0) find(params.xtime == 0)], [0 size(zmat(idx_exp,:,iSplit),1)+0.5],'k:','LineWidth',3);
        
        if ismember(iSplit,[7 8 9])
            set(gca, 'XTick', find(ismember(params.xtime,params.t_pre:0.5e6:params.t_post)), 'XTickLabels', (params.t_pre:0.5e6:params.t_post)/1e6,'FontSize', 20)
            xlabel('Time (s)','FontSize', 20)
        else
            set(gca, 'XTick', 1:1000:length(params.xtime), 'XTickLabels', [],'FontSize', 20)
        end
        
        if ismember(iSplit,[1 4 7])
            ylabel('Neuron','FontSize', 20)
            set(gca, 'YTick', [1 sum(idx_exp)], 'YTickLabels', [1 sum(idx_exp)],'FontSize', 20)
        else 
            set(gca, 'YTick', [], 'YTickLabels', [],'FontSize', 20)
        end
        set(gca,'box','off')
        title(params.triallabels{iSplit},'FontSize',15)
        xlim([find(params.xtime == -1e6) find(params.xtime == 1.5e6)]);
%         xlim([find(params.xtime == -2e6) find(params.xtime == 2e6)]);
    end
%     saveas(gcf,fullfile(params.savedir,'Zmean_9conditions.eps'))
%     saveas(gcf,fullfile(params.savedir,sprintf('Zmat_3x3_%s.eps',params.ExperimentLabels{iExp})))
%     export_fig(fullfile(params.savedir,sprintf('Zmat_3x3_%s.eps',params.ExperimentLabels{iExp})),'-painters')
%     export_fig(fullfile(params.savedir,sprintf('Zmat_3x3_%s.eps',params.ExperimentLabels{iExp})))
end

%make separate colorbar figure:
figure; set(gcf,'units','normalized','Position',[0.5 0.3 0.13 0.2],'color','w')
imagesc(zmat(idx_exp,:,iSplit),params.cscale); hold on;
c = colorbar;
c.Label.String = 'Z-scored firing rate';
cla

%% Make figure of the mean:
params.nExperiments = length(params.Experiments);
figure; set(gcf,'units','normalized','Position',[0.05 0.4 0.91 0.17],'color','w')
handles = cell(9,1);

for iExp = 1:params.nExperiments
    idx_exp = ismember(spikeData.session_ID,sessionData.session_ID(strcmp(sessionData.Experiment,params.Experiments(iExp))));
    
    for iSplit = 1:length(splits)
        iSub = (iExp-1)*3+mod(iSplit-1,3)+1;
        subplot(1,9,iSub)
        
        meantoplot      = nanmean(zmat(idx_exp,:,iSplit),1);
        errortoplot     = nanstd(zmat(idx_exp,:,iSplit),1)/sqrt(size(zmat(idx_exp,:,iSplit),1));
        h = shadedErrorBar(params.xtime,meantoplot,errortoplot,{params.triallines{iSplit},'markerfacecolor',params.trialcolors{iSplit},'LineWidth',3},0);
        if ~all(isnan(meantoplot))
            h.mainLine.Color = params.trialcolors{iSplit};   h.patch.FaceColor = params.trialcolors{iSplit};
            delete(h.edge(1)); delete(h.edge(2));
            handles{iSub}(end+1) = h.mainLine; hold all;
        end
    end
end

%Figure makeup:
for iSub = 1:9
    subplot(1,9,iSub)
    xlim([-1e6 1.5e6]);
    ylim([-0.3 1])
    set(gca, 'XTick', [-1e6 0 1.5e6], 'XTickLabels', [-1 0 1.5],'YTick',[-0.3 0 1],'FontSize', 10)
    plot([0 0],[-0.3 1],'k:','LineWidth',2) %Plot reference line at t=0, first lick
    h = patch([-0.3e6 0.3e6 0.3e6 -0.3e6],[-1 -1 1 1],[0.5 0.5 0.5]);
    if iSub==1 %For the first one plot the labels:
        xlabel('Time relative to first lick (s)','FontSize', 12)
        ylabel('Z-scored firing rate','FontSize', 12)
    end
    
    if iSub<4 %For first three plot the legend
        legend(handles{iSub},params.triallabels([1 4 7] + (iSub-1)),'FontSize',11,'Location','NorthWest'); legend boxoff
    end
    
end

%     export_fig(fullfile(params.savedir,'Zmean_9conditions.eps'))
% saveas(gcf,fullfile(params.savedir,'Zmean_9conditions.eps'))

%% Reshape data and average over timewindow around lick:
params.nSplits = 9;

timeidx     = params.xtime>-0.3e6 & params.xtime<0.3e6;

respmat     = NaN(params.nExperiments,params.nSplits,1000);

for iExp = 1:params.nExperiments
    idx_exp = ismember(spikeData.session_ID,sessionData.session_ID(strcmp(sessionData.Experiment,params.Experiments(iExp))));
    for iSplit = 1:params.nSplits
%         datatoplot(iExp,idx_reorder(iSplit),1:sum(idx_exp)) = nanmean(zmat(idx_exp,timeidx,iSplit),2); 
        respmat(iExp,iSplit,1:sum(idx_exp)) = nanmean(zmat(idx_exp,timeidx,iSplit),2); 
%         respmat(iExp,iSplit,idx_exp) = nanmean(zmat(idx_exp,timeidx,iSplit),2); 
    end
end


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

%% Figure 1: all responding conditions vs no-responding conditions split per cohort:

params.nSplits = 2;

idx_nolick      = [3 6 9];
idx_lick        = [1 2 4 5 7 8];
datatoplot = [];
datatoplot(:,1,:) = nanmean(respmat(:,idx_nolick,:),2);
datatoplot(:,2,:) = nanmean(respmat(:,idx_lick,:),2);
%Overwrite UST lick data because no licks to auditory lick spout:
idx_lickUST        = [2 5 8];
datatoplot(2,2,:) = nanmean(respmat(2,idx_lickUST,:),2);

meantoplot          = nanmean(datatoplot,3);
errortoplot         = nan(size(meantoplot));
for iExp = 1:params.nExperiments %For each experiment copute the mean and error by averaging expl var over neurons
    idx_exp                         = ismember(spikeData.session_ID,sessionData.session_ID(strcmp(sessionData.Experiment,params.Experiments(iExp))));
    errortoplot(iExp,:)             = nanstd(datatoplot(iExp,:,:),[],3) / sqrt(sum(idx_exp));
end

%% Make figure;
figure; hold all; set(gcf,'units','normalized','Position',[0.1141    0.3537    0.3349    0.4250],'color','w');

% violinplot(datatoplot)
h = violinplot(reshape(permute(datatoplot,[2 1 3]),6,1000)',repmat({'No Lick' 'Lick'},1,3),'ShowData',false,'ViolinAlpha',1);

for iExp = 1:params.nExperiments %For each experiment copute the mean and error by averaging expl var over neurons
    h(1,iExp*2-1).ViolinColor = params.colors_experiments{iExp};
    h(1,iExp*2).ViolinColor = params.colors_experiments{iExp};
end
set(gca,'YTick',-0.5:0.5:2)
xlim([0.5 6.5])
ylim([-0.5 2])

% %% OLD figure;
% figure; hold all; set(gcf,'units','normalized','Position',[0.1141    0.3537    0.3349    0.4250],'color','w');
% xpos = [0.8 1.2];
% for i = 1:size(meantoplot,1) %Give bars colors:
%     h = bar(xpos+i-1,meantoplot(i,:),0.8); %Plot bars with width 0.8
%     h.FaceColor = params.colors_experiments{i};
% end
% 
% % Calculating the width for each bar group
% for i = 1:size(meantoplot,1)
%     errorbar(xpos+i-1, meantoplot(i,:), errortoplot(i,:), 'k.', 'LineWidth',5); %show errorbar
% end
% set(gca,'XTick',[],'YTick',[0 0.3 0.6])

%% Statistical testing:

idx_nolick      = [3 6 9];
idx_lick        = [1 2 4 5 7 8];

Y_resp      = [squeeze(nanmean(nanmean(zmat(:,timeidx,idx_nolick),2),3)); squeeze(nanmean(nanmean(zmat(:,timeidx,idx_lick),2),3))];

X_lick      = [ones(nNeurons,1); ones(nNeurons,1)*2];

tbl         = table(Y_resp,X_lick,[X_coh; X_coh],[G_mou; G_mou],'VariableNames',{'Activity','Lick','Cohort','Mouse'}); %Create table for mixed model

lme         = fitlme(tbl,'Activity~Lick+Cohort+Lick*Cohort'); %construct linear mixed effects model with fixed effect of temporal window and random intercept for different mice
stats       = dataset2table(anova(lme,'DFMethod','Satterthwaite')); %Perform ANOVA on model and output as matrix

fprintf('Fixed effect of Licking: ')
fprintf('F(%d,%2.0f)=%1.2f, p=%1.2e; \n',stats{2,3},stats{2,4},stats{2,2},stats{2,5})
fprintf('Fixed effect of Cohort: ')
fprintf('F(%d,%2.0f)=%1.2f, p=%1.3f; \n',stats{3,3},stats{3,4},stats{3,2},stats{3,5})
fprintf('Interaction effect Licking * Cohort: ')
fprintf('F(%d,%2.0f)=%1.2f, p=%1.3e; \n',stats{4,3},stats{4,4},stats{4,2},stats{4,5})

contrast = [0 0 0 0 1 0]; %compare interaction term NE with MST/UST
fprintf('Posthoc comparison: lick modulation NE vs UST/MST: \n')
[p,F,DF1,DF2] = coefTest(lme,contrast);
fprintf('F(%d,%d)=%1.2f, p=%1.3e\n',DF1,DF2,F,p)

tempfile = fullfile('SourceData_FigS3i_Activity_Lick_Cohort.xlsx');
writetable(tbl,tempfile)

% p = stats{2,5};
% if p<0.05
%     sigstar([0.5 1],p)
% end

% OLD Statistical testing:
% groups              = repmat(1:params.nExperiments*params.nSplits,1000,1); groups = reshape(groups,params.nExperiments*params.nSplits*1000,1);
% datatotest          = reshape(permute(datatoplot,[3 2 1]),params.nExperiments*params.nSplits*1000,1); %reshape to one column vector
% groups              = groups(~isnan(datatotest)); %filter out nans
% datatotest          = datatotest(~isnan(datatotest)); %filter out nans
% %perform kruskal wallis nonparametric anova:
% [p,table,stats]     = kruskalwallis(datatotest,groups,'off');
% comptable           = multcompare(stats,'display','off','alpha',params.alpha,'ctype',params.posthoctest); %do posthoc
% 
% comptable           = comptable(comptable(:,end)<params.alpha,:); %Filter only significant
% xpos = [xpos xpos+1 xpos+2];
% sigstar(mat2cell(xpos(comptable(:,1:2)),ones(size(comptable,1),1)),comptable(:,end)) %use sigstar function to identify
% set(gca,'XTick',[],'YTick',[0 0.25 0.5 0.75])
% ylim([-0.2 inf])

%%

%% Figure 2: all hits vs all errors (within responding trials)

params.nSplits  = 6;

datatoplot      = NaN(6,1000);

datatoplot(1,:) = nanmean(respmat(2,[2 8],:),2); %visual errors (FA and incorrect) for UST
datatoplot(2,:) = nanmean(respmat(2,5,:),2); %visual hits UST, etc...

datatoplot(3,:) = nanmean(respmat(3,[2 8],:),2);
datatoplot(4,:) = nanmean(respmat(3,5,:),2);

datatoplot(5,:) = nanmean(respmat(3,[4 7],:),2);
datatoplot(6,:) = nanmean(respmat(3,1,:),2);

meantoplot      = nanmean(datatoplot,2);
errortoplot     = nan(size(meantoplot));
for iSplit = 1:params.nSplits %For each experiment copute the mean and error by averaging expl var over neurons
    errortoplot(iSplit,:)             = nanstd(datatoplot(iSplit,:),[],2) / sqrt(sum(~isnan(datatoplot(iSplit,:))));
end

% %% Make old figure;
% figure; hold all; set(gcf,'units','normalized','Position',[0.1141    0.3537    0.3349    0.4250],'color','w');
% xpos = [0.8 1.2];
% h = bar(xpos,meantoplot([1 2],:),0.8); %Plot bars with width 0.8
% h.FaceColor = params.colors_modalities{2};
% errorbar(xpos, meantoplot([1 2],:), errortoplot([1 2],:), 'k.', 'LineWidth',5); %show errorbar
% 
% h = bar(xpos+1,meantoplot([3 4],:),0.8); %Plot bars with width 0.8
% h.FaceColor = params.colors_modalities{2};
% errorbar(xpos+1, meantoplot([3 4],:), errortoplot([3 4],:), 'k.', 'LineWidth',5); %show errorbar
% 
% h = bar(xpos+2,meantoplot([5 6],:),0.8); %Plot bars with width 0.8
% h.FaceColor = params.colors_modalities{1};
% errorbar(xpos+2, meantoplot([5 6],:), errortoplot([5 6],:), 'k.', 'LineWidth',5); %show errorbar
% 
% set(gca,'XTick',[],'YTick',[0 0.25 0.5 0.75])
% ylim([0 0.75])

%% Make figure;
figure; hold all; set(gcf,'units','normalized','Position',[0.1141    0.3537    0.3349    0.4250],'color','w');

% violinplot(datatoplot)
h = violinplot(datatoplot',repmat({'*Incorr' 'Corr'},1,3),'ShowData',false,'ViolinAlpha',1,'Bandwidth',0.1);

% for iExp = 1:params.nExperiments %For each experiment copute the mean and error by averaging expl var over neurons
h(1,1).ViolinColor = params.colors_modalities{2};
h(1,2).ViolinColor = params.colors_modalities{2};
h(1,3).ViolinColor = params.colors_modalities{2};
h(1,4).ViolinColor = params.colors_modalities{2};
h(1,5).ViolinColor = params.colors_modalities{1};
h(1,6).ViolinColor = params.colors_modalities{1};

%Statistics:
iExp            = 2;
idx_exp         = ismember(spikeData.session_ID,sessionData.session_ID(strcmp(sessionData.Experiment,params.Experiments(iExp))));
Y_resp          = [squeeze(nanmean(nanmean(zmat(idx_exp,timeidx,[2 8]),2),3)); squeeze(nanmean(nanmean(zmat(idx_exp,timeidx,5),2),3))];
X_corr          = [ones(sum(idx_exp),1); ones(sum(idx_exp),1)*2];
tbl             = table(Y_resp,X_corr,[G_mou(idx_exp); G_mou(idx_exp)],'VariableNames',{'Activity','Correct','Mouse'}); %Create table for mixed model
lme             = fitlme(tbl,'Activity~Correct+(1|Mouse)'); %construct linear mixed effects model with fixed effect of temporal window and random intercept for different mice
stats           = dataset2table(anova(lme,'DFMethod','Satterthwaite')); %Perform ANOVA on model and output as matrix
fprintf('UST - visual lick spout: \n')
fprintf('F(%d,%2.0f) = %1.2f, p=%1.3f; \n',stats{2,3},stats{2,4},stats{2,2},stats{2,5})
fprintf('(Linear Mixed Model ANOVA)\n')
sigstar(xpos,stats{2,5});

tempfile = fullfile('SourceData_FigS3j_Activity_HitError_UST_Vis.xlsx');
writetable(tbl,tempfile)

iExp            = 3;
idx_exp         = ismember(spikeData.session_ID,sessionData.session_ID(strcmp(sessionData.Experiment,params.Experiments(iExp))));
Y_resp          = [squeeze(nanmean(nanmean(zmat(idx_exp,timeidx,[2 8]),2),3)); squeeze(nanmean(nanmean(zmat(idx_exp,timeidx,5),2),3))];
X_corr          = [ones(sum(idx_exp),1); ones(sum(idx_exp),1)*2];
tbl             = table(Y_resp,X_corr,[G_mou(idx_exp); G_mou(idx_exp)],'VariableNames',{'Activity','Correct','Mouse'}); %Create table for mixed model
lme             = fitlme(tbl,'Activity~Correct+(1|Mouse)'); %construct linear mixed effects model with fixed effect of temporal window and random intercept for different mice
stats           = dataset2table(anova(lme,'DFMethod','Satterthwaite')); %Perform ANOVA on model and output as matrix
fprintf('MST - visual lick spout: \n')
fprintf('F(%d,%2.0f) = %1.2f, p=%1.2e; \n',stats{2,3},stats{2,4},stats{2,2},stats{2,5})
fprintf('(Linear Mixed Model ANOVA)\n')
sigstar(xpos+2,stats{2,5});

tempfile = fullfile('SourceData_FigS3j_Activity_HitError_MST_Vis.xlsx');
writetable(tbl,tempfile)

iExp            = 3;
idx_exp         = ismember(spikeData.session_ID,sessionData.session_ID(strcmp(sessionData.Experiment,params.Experiments(iExp))));
Y_resp          = [squeeze(nanmean(nanmean(zmat(idx_exp,timeidx,[4 7]),2),3)); squeeze(nanmean(nanmean(zmat(idx_exp,timeidx,1),2),3))];
X_corr          = [ones(sum(idx_exp),1); ones(sum(idx_exp),1)*2];
tbl             = table(Y_resp,X_corr,[G_mou(idx_exp); G_mou(idx_exp)],'VariableNames',{'Activity','Correct','Mouse'}); %Create table for mixed model
lme             = fitlme(tbl,'Activity~Correct+(1|Mouse)'); %construct linear mixed effects model with fixed effect of temporal window and random intercept for different mice
stats           = dataset2table(anova(lme,'DFMethod','Satterthwaite')); %Perform ANOVA on model and output as matrix
fprintf('MST - audio lick spout: \n')
fprintf('F(%d,%2.0f) = %1.2f, p=%1.3f; \n',stats{2,3},stats{2,4},stats{2,2},stats{2,5})
fprintf('(Linear Mixed Model ANOVA)\n')
sigstar(xpos+4,stats{2,5});

tempfile = fullfile('SourceData_FigS3j_Activity_HitError_MST_Aud.xlsx');
writetable(tbl,tempfile)

set(gca,'YTick',-0.5:0.5:2)
xlim([0.5 6.5])
ylim([-1 2.5])
xpos = [1 2];

% fprintf('Order of testing: Visual spout (UST), Visual Spout (MST), Auditory spout (MST):\n')
% for i = 1:3
%     idx = i*2-1 : i*2;
%     p = signrank(datatoplot(idx(1),:),datatoplot(idx(2),:));
%     p = p*3;
%     sigstar(xpos+i-1,p);
%     fprintf('\nWilcoxon signed rank test Errors vs Hits - %d neurons: %1.2e\n',sum(~isnan(datatoplot(idx(2),:))),p)
% end


