%% This script analyzes lick related activity
% This script aligns V1 neural activity to stimulus onset and compares early and late activity dynamics
% Reproduces panels from main figure 2
% Oude Lohuis et al. 2022 Nat Comms

%% Reset all
startover

%% Parameter settings:
params                      = params_histresponse_coding; %Parameters for PSTH (All time is in microseconds)
params.conv_win             = 'chg';

params.Experiments          = {'ChangeDetectionConflictDecor' 'VisOnlyTwolevels' 'ChangeDetectionConflict' }; %Which versions of the task to load data from
params.ExperimentLabels     = {'NE' 'UST' 'MST'}; %Labels for each task version
params.nExperiments         = length(params.Experiments);

params.AlignOn              = 'stimChange';      %On which timestamp to align as t=0

params.clipzscore           = 25;

%parameters for display of zscore matrix:
params.colormap             = 'parula';
params.cscale               = [-0.5 2.2];

params.minTrialCond         = 3;
% params.cscale               = [-1 2];

params.area                 = 'V1'; %Filter only V1 data

params.twin_resp_start      = 0e6;
params.twin_resp_stop       = 0.5e6;

params                      = MOL_getColors_CHDET(params);

params.trialcategories      = 'ProbeVis';

params.labels_ztrials       = {'MISS' 'HIT' 'CR' 'FA'};
params.lines_ztrials        = {'-' '-' '-' '-'};

params.lines_experiments    = {'-' '-' '-'};

params.savedir              = 'E:\Documents\PhD\Figures\Project CHDET\Results - 2nd Bump\2 Neural\Zmat_AlignStim';

%% Load dataset:
load('E:\Matlab\MOL_Analysis\AOudeLohuisetal_2022_NatComms\2Neural\Data2_1.mat')
fprintf('Dataset: %d sessions, %d trials, %d neurons\n',length(sessionData.session_ID),length(trialData.session_ID),length(spikeData.session_ID));

%% Show raster plot of example neurons showing multiplexing of sensory and decision components in single V1 neurons
cell_IDs                = {};
%Naive examples:
cell_IDs{end+1}         = '10092019030831134';

%Alternatives:
% cell_IDs{end+1}     = '10082019030831 218';
% cell_IDs{end+1}     = '10092019030721318';

% UST example:
cell_IDs{end+1}         = '20282019121311450'; %UST

%Alternatives:
% cell_IDs{end+1}     = '20292019121211226'; %UST
% cell_IDs{end+1}     = '20292019121211218'; %UST
% cell_IDs{end+1}     = '20342019121931076'; %UST
% cell_IDs{end+1}     = '20342019121931054'; %UST
% cell_IDs{end+1}     = '20292019121211213'; %UST
% cell_IDs{end+1}     = '20282019121311482'; %UST
% cell_IDs{end+1}     = '20352019121931087'; %UST

% MST Example:
cell_IDs{end+1}         = '20122018081431146'; %MST

%Alternatives:
% cell_IDs{end+1}     = '20032018020821278'; %MST
% cell_IDs{end+1}     = '20122018081531198';
% cell_IDs{end+1}     = '20122018081531204';
% cell_IDs{end+1}     = '20222019062811204';

cell_IDs_zlabels        = {'a' 'b' 'c'};

outputmat_sign          = [];
params.AUC_varselec     = [];

params.export_fig       = 1;

MOL_plotRaster(sessionData,trialData,spikeData,cell_IDs,params)

%% Stratify reaction times between FA and visual hits:

figure; set(gcf,'units','normalized','Position',[0.1 0.6 0.48 0.25],'color','w')

splits{1}               = ismember(trialData.trialType,{'X'}) & trialData.visualOriChangeNorm==3 & trialData.vecResponse==2;
splits{2}               = ismember(trialData.trialType,{'P'}) & trialData.vecResponse==2;

iExp = 2;
idx_exp = ismember(trialData.session_ID,sessionData.session_ID(strcmp(sessionData.Experiment,params.Experiments(iExp))));
binedges = 0:25e3:1000e3;
nBins = numel(binedges)-1;

subplot(1,2,1); hold all;
title('UST')
RThist_UST_HIT = histcounts(trialData.responseLatency(idx_exp & splits{1}),binedges,'normalization','pdf');
plot(binedges(1:end-1)+12.5e3,RThist_UST_HIT,'-','Color',params.colors_experiments{2})
RThist_UST_FA = histcounts(trialData.responseLatency(idx_exp & splits{2}),binedges,'normalization','pdf');
plot(binedges(1:end-1)+12.5e3,RThist_UST_FA,'k:','Color',params.colors_experiments{2})

iExp = 3;
idx_exp = ismember(trialData.session_ID,sessionData.session_ID(strcmp(sessionData.Experiment,params.Experiments(iExp))));
subplot(1,2,2); hold all;
title('MST')
binedges = 0:25e3:1000e3;
RThist = histcounts(trialData.responseLatency(idx_exp & splits{1}),binedges,'normalization','pdf');
plot(binedges(1:end-1)+12.5e3,RThist,'-','Color',params.colors_experiments{3})
RThist = histcounts(trialData.responseLatency(idx_exp & splits{2}),binedges,'normalization','pdf');
plot(binedges(1:end-1)+12.5e3,RThist,'k:','Color',params.colors_experiments{3})

iExp = 2;
idx_exp = ismember(trialData.session_ID,sessionData.session_ID(strcmp(sessionData.Experiment,params.Experiments(iExp))));
frac = RThist_UST_FA ./ RThist_UST_HIT;
idx_all_rem = [];
for iBin = 1:nBins
    if frac(iBin)>1
       idx_RT       = trialData.responseLatency>binedges(iBin) & trialData.responseLatency<=binedges(iBin+1);
       idx          = find(splits{2} & idx_exp & idx_RT);
%        idx          = splits{1} & idx_exp & idx_RT;
       nSubselec    = floor(numel(idx) * 1/frac(iBin));
       
       subselec     = idx(randperm(numel(idx),nSubselec));
       idx_all_rem  = [idx_all_rem; idx(~ismember(idx,subselec))]; %#ok<AGROW>
    end
end

trialData.trialType(idx_all_rem) = {'R'};

%% Main loop to get psth matrix:
params.nSplits          = 4;
params.zscore           = 1;

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
    splits{1}               = ismember(temptrialData.trialType,{'X'}) & temptrialData.visualOriChangeNorm==3 & temptrialData.vecResponse==3;
    splits{2}               = ismember(temptrialData.trialType,{'X'}) & temptrialData.visualOriChangeNorm==3 & temptrialData.vecResponse==2;
    splits{3}               = ismember(temptrialData.trialType,{'P'}) & temptrialData.vecResponse==3;
    splits{4}               = ismember(temptrialData.trialType,{'P'}) & temptrialData.vecResponse==2;
    
    for iSplit = 1:params.nSplits %Store the mean response for each of these splits
        if sum(splits{iSplit})>=params.minTrialCond
            zmat(iNeuron,:,iSplit) = mean(hist_mat(splits{iSplit},:),1);
        end
    end
end

%% Sort zmat
% meanresp                = mean(zmat(:,params.xtime>params.twin_resp_start & params.xtime<=params.twin_resp_stop,3),2) + mean(zmat(:,params.xtime>params.twin_resp_start & params.xtime<=params.twin_resp_stop,4),2);
meanresp                = mean(zmat(:,params.xtime>params.twin_resp_start & params.xtime<=params.twin_resp_stop,1),2) + mean(zmat(:,params.xtime>params.twin_resp_start & params.xtime<=params.twin_resp_stop,2),2);

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
    title(params.ExperimentLabels{iExp})
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
                yt = get(gca,'YTick'); yt(end+1) = sum(idx_all)-temp; %#ok<SAGROW>
                ytl = get(gca,'YTickLabels'); ytl(end+1) = cell_IDs_zlabels(i); %#ok<SAGROW>
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

%% Make figure of the mean:
for iExp = 1:params.nExperiments
    idx_exp = ismember(spikeData.session_ID,sessionData.session_ID(strcmp(sessionData.Experiment,params.Experiments(iExp))));
    idx_all = idx_exp & ~idx_nan;
    
    figure; set(gcf,'units','normalized','Position',[0.4 0.3 0.32 0.38],'color','w')
    title(params.ExperimentLabels{iExp})
    
    %     handles = NaN(params.nSplits,1);
    handles = [];
    for iSplit = 1:params.nSplits
        meantoplot = nanmean(zmat(idx_all,:,iSplit),1);
        errortoplot = nanstd(zmat(idx_all,:,iSplit),1)/sqrt(size(zmat(idx_all,:,iSplit),1));
        plot([0 0], [-2 5],'k:','LineWidth',3);

        if ~all(isnan(meantoplot))
            h = shadedErrorBar(params.xtime,meantoplot,errortoplot,{'-k','markerfacecolor',params.colors_ztrials{iSplit},'LineWidth',3},1);
            h.mainLine.Color = params.colors_ztrials{iSplit};    h.patch.FaceColor = params.colors_ztrials{iSplit};
            delete(h.edge(1)); delete(h.edge(2)); 
            handles(end+1) = h.mainLine; hold all; %#ok<SAGROW>
        end
    end
%     set(gca, 'XTick', find(ismember(params.xtime,-3e6:0.5e6:3e6)), 'XTickLabels', params.xtime(ismember(params.xtime,(-3e6:0.5e6:3e6)))/1e6,'FontSize', 20)
    set(gca, 'XTick', -3e6:0.5e6:3e6, 'XTickLabels',(-3e6:0.5e6:3e6)/1e6,'FontSize', 20)
    set(gca, 'YTick', [0 0.5 1 1.5], 'FontSize', 20)
    xlim([-0.5e6 1.5e6]);
    ylim([-0.2 1.5])
    
    ylabel('Z-scored firing rate','FontSize', 20)
    xlabel('Time (s)','FontSize', 20)
    legend(handles,params.labels_ztrials); legend boxoff
end


%%
VisGranularUpper        = 400; %Upper boundary of layer IV in V1
VisGranularLower        = 550; %Lower boundary of layer IV in V1

spikeData.Layer((spikeData.ChannelY)<=VisGranularUpper,1) = deal({'SG'});
spikeData.Layer((spikeData.ChannelY)>VisGranularUpper & (spikeData.ChannelY)<VisGranularLower,1) = deal({'G'});
spikeData.Layer((spikeData.ChannelY)>=VisGranularLower,1) = deal({'IG'});

spikeData.Layer((spikeData.ChannelY)<0,1) = deal({'YY'});
spikeData.Layer((spikeData.ChannelY)>1000,1) = deal({'XX'});

%% Make figure of the mean:
iExp = 3;
figure; set(gcf,'units','normalized','Position',[0.06 0.3 0.64 0.25],'color','w')

idx_exp = ismember(spikeData.session_ID,sessionData.session_ID(ismember(sessionData.Experiment,params.Experiments([2 3]))));

datatotest = NaN(1000,2,3); %neurons x early/late x layer

layers = {'SG' 'G' 'IG'};
for iLayer = 1:3
    idx_layer = strcmp(spikeData.Layer,layers{iLayer});
    idx_all = idx_exp & ~idx_nan & idx_layer;
    
    subplot(1,3,iLayer);

    handles     = [];
    
    idx_early   = params.xtime>0 & params.xtime<=200e3;
    idx_late    = params.xtime>200e3 & params.xtime<=1000e3;
    
    datatotest(1:sum(idx_all),1,iLayer) = squeeze(nanmax(zmat(idx_all,idx_early,2),[],2));
    datatotest(1:sum(idx_all),2,iLayer) = squeeze(nanmax(zmat(idx_all,idx_late,2),[],2));
    
    for iSplit = [1 2]
        meantoplot = nanmean(zmat(idx_all,:,iSplit),1);
        errortoplot = nanstd(zmat(idx_all,:,iSplit),1)/sqrt(size(zmat(idx_all,:,iSplit),1));
        plot([0 0], [-2 5],'k:','LineWidth',3);

        if ~all(isnan(meantoplot))
            h = shadedErrorBar(params.xtime,meantoplot,errortoplot,{'-k','markerfacecolor',params.colors_ztrials{iSplit},'LineWidth',3},0);
            h.mainLine.Color = params.colors_ztrials{iSplit};    h.patch.FaceColor = params.colors_ztrials{iSplit};
            delete(h.edge(1)); delete(h.edge(2)); %h.edge(1).Color = [1 1 1];  h.edge(2).Color =  [1 1 1]; h.edge(1).LineWidth = 0.1;  h.edge(2).LineWidth = 0.1;
            handles(end+1) = h.mainLine; hold all; %#ok<SAGROW>
        end
    end
%     set(gca, 'XTick', find(ismember(params.xtime,-3e6:0.5e6:3e6)), 'XTickLabels', params.xtime(ismember(params.xtime,(-3e6:0.5e6:3e6)))/1e6,'FontSize', 20)
    set(gca, 'XTick', -3e6:0.5e6:3e6, 'XTickLabels',(-3e6:0.5e6:3e6)/1e6,'FontSize', 20)
    set(gca, 'YTick', [0 1 2], 'FontSize', 20)
    xlim([-0.5e6 1.5e6]);
    ylim([-0.3 3])
    
    title(layers{iLayer})
    ylabel('Z-scored firing rate','FontSize', 20)
    xlabel('Time (s)','FontSize', 20)
    legend(handles,params.labels_ztrials(1:2)); legend boxoff
end

%% Multi-level statistics: 

X_coh           = NaN(nNeurons,1);
X_coh(ismember(spikeData.session_ID,sessionData.session_ID(ismember(sessionData.Experiment,params.Experiments{1}))))         = 1;
X_coh(ismember(spikeData.session_ID,sessionData.session_ID(ismember(sessionData.Experiment,params.Experiments{2}))))         = 2;
X_coh(ismember(spikeData.session_ID,sessionData.session_ID(ismember(sessionData.Experiment,params.Experiments{3}))))         = 3;

idx_early       = params.xtime>0 & params.xtime<=200e3;
idx_late        = params.xtime>200e3 & params.xtime<=1000e3;

G_mou           = NaN(nNeurons,1);
uMice           = unique(sessionData.mousename);
for iMouse = 1:length(uMice)
    G_mou(ismember(spikeData.session_ID,sessionData.session_ID(strcmp(sessionData.mousename,uMice{iMouse})))) = iMouse;
end

%% Comparing early and late for each layer:

meantoplot = squeeze(nanmean(datatotest,1));

errortoplot = squeeze(nanstd(datatotest,[],1));
errortoplot = errortoplot ./ repmat(sqrt(squeeze(sum(~isnan(datatotest(:,1,:)),1)))',2,1);

sum(~isnan(datatotest(:,1,iLayer)))

Y_early         = squeeze(nanmean(zmat(:,idx_early,2),2));
Y_late          = squeeze(nanmean(zmat(:,idx_late,2),2));

figure; set(gcf,'units','normalized','Position',[0.26 0.45 0.16 0.18],'color','w'); hold all;

G_mou_2 = [G_mou; G_mou];
X_coh_2 = [X_coh; X_coh];

for iLayer = 1:3
    subplot(1,3,iLayer)
%     errorbar([0.5 1]+iLayer-1,meantoplot(:,iLayer),errortoplot(:,iLayer),'o-','Color',colors(iLayer,:),'MarkerSize',10,'MarkerFaceColor',[0 0 0],'LineWidth',2)
%     errorbar([0.5 1],meantoplot(:,iLayer),errortoplot(:,iLayer),'o-','Color',[0 0 0],'MarkerSize',8,'MarkerFaceColor',colors(iLayer,:),'LineWidth',1.5)
    errorbar([0.5 1],meantoplot(:,iLayer),errortoplot(:,iLayer),'o-','Color',[0 0 0],'MarkerSize',8,'MarkerFaceColor',[0.6 0.6 0.6],'LineWidth',1.5)
    Y               = [Y_early; Y_late];
    X_Layer_2       = [strcmp(spikeData.Layer,layers{iLayer}); strcmp(spikeData.Layer,layers{iLayer})];
    X_window        = [ones(size(Y_early)); ones(size(Y_late))*2];
    idx             = ismember(X_coh_2,[2 3]) & X_Layer_2==1; %only compare UST and MST
    tbl             = table(Y(idx),X_window(idx),G_mou_2(idx),'VariableNames',{'Activity','Window','Mouse'}); %Create table for mixed model
    
    lme             = fitlme(tbl,'Activity~Window+(1|Mouse)'); %construct linear mixed effects model with fixed effect of temporal window and random intercept for different mice
    stats           = dataset2table(anova(lme,'DFMethod','Satterthwaite')); %Perform ANOVA on model and output as matrix
    fprintf('%s: ',layers{iLayer})
    fprintf('F(%d,%2.0f)=%1.2f, p=%1.2f; \n',stats{2,3},stats{2,4},stats{2,2},stats{2,5})
    p = stats{2,5};
    if p<0.05
        sigstar([0.5 1],p)
    end
    ylim([1.3 4.5])
    %     ylim([-0.5 2])
    xlim([0.25 1.25])
    set(gca,'YTick',[1.5 4],'XTick',[0.5 1],'XTickLabel',{'Early' 'Late'},'XTickLabelRotation',45)
    
    tempfile = fullfile(sprintf('SourceData_Fig2f_Early_Late_%s.xlsx',layers{iLayer}));
    writetable(tbl,tempfile)
end
fprintf('Linear Mixed Model ANOVA\n')



%% Comparing late hit/miss modulation per layer:

Y_late_mod      = squeeze(nanmean(zmat(:,idx_late,2),2)) - squeeze(nanmean(zmat(:,idx_late,1),2));

X_Layer         = spikeData.Layer;

idx             = ismember(X_coh,[2 3]) & ismember(spikeData.Layer,layers); %index of only UST and MST

tbl             = table(Y_late_mod(idx),X_Layer(idx),G_mou(idx),'VariableNames',{'Modulation','Layer','Mouse'}); %Create table for mixed model
lme             = fitlme(tbl,'Modulation~Layer+(1|Mouse)'); %construct linear mixed effects model with fixed effect random intercept for different mice
stats           = dataset2table(anova(lme,'DFMethod','Satterthwaite')); %Perform ANOVA on model and output as matrix
fprintf('Late modulation, fixed effect of layer: (Linear Mixed Model)\n')
fprintf('(F(%d,%2.0f) = %1.2f, p=%1.3f, ANOVA)\n',stats{2,3},stats{2,4},stats{2,2},stats{2,5})

contrasts = {[1 -1 0] [0 1 -1] [1 0 -1]};
contrastlabels = {'IG vs G' 'G vs SG' 'IG vs SG'};

fprintf('Posthoc comparison:\n')
for iC = 1:3
    [p,F,DF1,DF2] = coefTest(lme,contrasts{iC});
    fprintf('(%s):F(%d,%d)=%1.2f, p=%1.2f\n',contrastlabels{iC},DF1,DF2,F,p)
end

tempfile = fullfile(sprintf('SourceData_Fig2f_HitMiss_DiffAcrossLayers.xlsx'));
writetable(tbl,tempfile)

%% Comparing early activity for hits and misses for each layer:
%(For reviewer)
Y_hits         = squeeze(nanmean(zmat(:,idx_early,2),2));
Y_miss         = squeeze(nanmean(zmat(:,idx_early,1),2));

G_mou_2 = [G_mou; G_mou];
X_coh_2 = [X_coh; X_coh];

%overall difference:
Y               = [Y_hits; Y_miss];
X_window        = [ones(size(Y_hits)); ones(size(Y_miss))*2];
idx             = ismember(X_coh_2,[2 3]); %only compare UST and MST
tbl             = table(Y(idx),X_window(idx),G_mou_2(idx),'VariableNames',{'Activity','Hits','Mouse'}); %Create table for mixed model

lme             = fitlme(tbl,'Activity~Hits+(1|Mouse)'); %construct linear mixed effects model with fixed effect of temporal window and random intercept for different mice
stats           = dataset2table(anova(lme,'DFMethod','Satterthwaite')); %Perform ANOVA on model and output as matrix
fprintf('Early vs late activity (all layers): (Linear Mixed Model)\n')
fprintf('(F(%d,%2.0f) = %1.2f, p=%1.2f, ANOVA)\n',stats{2,3},stats{2,4},stats{2,2},stats{2,5})
p = stats{2,5};
    
for iLayer = 1:3
    Y               = [Y_hits; Y_miss];
    X_Layer_2       = [strcmp(spikeData.Layer,layers{iLayer}); strcmp(spikeData.Layer,layers{iLayer})];
    X_window        = [ones(size(Y_hits)); ones(size(Y_miss))*2];
    idx             = ismember(X_coh_2,[2 3]) & X_Layer_2==1; %only compare UST and MST
    tbl             = table(Y(idx),X_window(idx),G_mou_2(idx),'VariableNames',{'Activity','Hits','Mouse'}); %Create table for mixed model
    
    lme             = fitlme(tbl,'Activity~Hits+(1|Mouse)'); %construct linear mixed effects model with fixed effect of temporal window and random intercept for different mice
    stats           = dataset2table(anova(lme,'DFMethod','Satterthwaite')); %Perform ANOVA on model and output as matrix
    fprintf('Early vs late activity (%s): (Linear Mixed Model)\n',layers{iLayer})
    fprintf('(F(%d,%2.0f) = %1.2f, p=%1.2f, ANOVA)\n',stats{2,3},stats{2,4},stats{2,2},stats{2,5})
    p = stats{2,5};
end

%% OLD Statistics
% meantoplot = squeeze(nanmean(datatotest,1));
% 
% errortoplot = squeeze(nanstd(datatotest,[],1));
% errortoplot = errortoplot ./ repmat(sqrt(squeeze(sum(~isnan(datatotest(:,1,:)),1)))',2,1);
% 
% sum(~isnan(datatotest(:,1,iLayer)))
% 
% figure; set(gcf,'units','normalized','Position',[0.26 0.45 0.16 0.18],'color','w'); hold all;
% 
% % colors = {[0.4 0.4 0.4] [0.4 0.4 0.4] [0.4 0.4 0.4]};
% colors = getPyPlot_cMap('BuPu',5);
% 
% for iLayer = 1:3
%     subplot(1,3,iLayer)
% %     errorbar([0.5 1]+iLayer-1,meantoplot(:,iLayer),errortoplot(:,iLayer),'o-','Color',colors(iLayer,:),'MarkerSize',10,'MarkerFaceColor',[0 0 0],'LineWidth',2)
% %     errorbar([0.5 1],meantoplot(:,iLayer),errortoplot(:,iLayer),'o-','Color',[0 0 0],'MarkerSize',8,'MarkerFaceColor',colors(iLayer,:),'LineWidth',1.5)
%     errorbar([0.5 1],meantoplot(:,iLayer),errortoplot(:,iLayer),'o-','Color',[0 0 0],'MarkerSize',8,'MarkerFaceColor',[0.6 0.6 0.6],'LineWidth',1.5)
% %     p = signrank(datatotest(:,1,iLayer),datatotest(:,2,iLayer));
% %     if p<0.05
% %         sigstar([0.5 1],p)
% %     end
% %     fprintf('Wilcoxon signed rank test early vs late response %s - %d neurons: %1.3f\n',layers{iLayer},sum(~isnan(datatotest(:,1,iLayer))),p)
%     ylim([1.3 4])
%     xlim([0.25 1.25])
%     set(gca,'YTick',[1.5 3],'XTick',0.75:2.75,'XTickLabel',{'SG' 'G' 'IG'})
%     set(gca,'YTick',[1.5 3],'XTick',[0.5 1],'XTickLabel',{'Early' 'Late'},'XTickLabelRotation',45)
% end
