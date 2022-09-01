%% This script analyzes noise correlations over time for different trial types
% Aligned to stimulus onset or to first lick
% Script reproduces panels from Figure 7
% Oude Lohuis et al. 2022 Nat Comms

%% Reset all
startover

%% Parameter settings:
params                  = params_histresponse(); %Parameters for PSTH (All time is in microseconds)

params.Experiments      = {'ChangeDetectionConflictDecor' 'VisOnlyTwolevels' 'ChangeDetectionConflict' }; %Which versions of the task to load data from
params.ExperimentLabels = {'NE' 'UST' 'MST'}; %Labels for the different experiments
params.nExperiments     = length(params.Experiments);

params.AlignOn          = 'stimChange';      %On which timestamp to align as t=0

params.minNtrialscond   = 10;

params.area             = 'V1';

params.savedir          = 'E:\Documents\PhD\Figures\Project CHDET\Results - 2nd Bump\6 Noisecorrelations';

params                  = MOL_getColors_CHDET(params); %Get colors and labels:

%% Load dataset:
load('E:\Matlab\MOL_Analysis\AOudeLohuisetal_2022_NatComms\2Neural\Data2_1.mat')
fprintf('Dataset: %d sessions, %d trials, %d neurons\n',length(sessionData.session_ID),length(trialData.session_ID),length(spikeData.session_ID));






%% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Analysis 1:
% Show noise correlations over time for each of the trial types aligned to stimulus onset
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%% Trial type splits:

%Parameters for trial splits:
params.labels_splits            = {'AuMiss' 'AuHit' 'VisMiss' 'VisHit' 'CR' 'FA'};
params.lines_splits             = {'-.' '-' '-.' '-' '-.' '-'};

% Get indices of the trials per trial type:
splits              = {};
splits{end+1}       = ismember(trialData.trialType,{'Y'}) & trialData.vecResponse==3;
splits{end+1}       = ismember(trialData.trialType,{'Y'}) & trialData.vecResponse==1;
splits{end+1}       = ismember(trialData.trialType,{'X'}) & trialData.vecResponse==3;
splits{end+1}       = ismember(trialData.trialType,{'X'}) & trialData.vecResponse==2;
splits{end+1}       = ismember(trialData.trialType,{'P'}) & trialData.vecResponse==3;
splits{end+1}       = ismember(trialData.trialType,{'P'}) & ismember(trialData.vecResponse,[1 2]);

%% Compute noise correlations: (this will take time)
%Set parameters:
params                      = params_histresponse_coding(params);
params.binsize              = 25e3;
params.smoothing            = 1;
params.conv_win             = 'chg';       %Type of window used for smoothing {flat, gaussian)
params.conv_sigma           = 0.100e6;           %sd of gaussian window for smoothing
params.conv_twin            = round((9*params.conv_sigma)/params.binsize)*params.binsize;         %Window size for smoothing

%Construct bin edges and time axis
params.edges                = [params.t_pre:params.binsize:params.t_post] - params.binsize/2;                    %#ok<NBRAK> Define edges of the bins
params.xtime                = params.t_pre:params.binsize:params.t_post-params.binsize;                        %Define time axis
params.nTimebins            = length(params.xtime); %number of time bins

params.AlignOn              = 'stimChange';

%Some measures of dimensionality:
nTotalNeurons               = length(spikeData.session_ID);
nSessions                   = length(sessionData.session_ID);
nConditions                 = numel(splits);
noisecorrmat                = NaN(nConditions,nTotalNeurons,nTotalNeurons,params.nTimebins); %init output matrix
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
    
    %subtract mean response for each condition:
    for iCond = 1:nConditions
        trialidx                = splits{iCond}(strcmp(trialData.session_ID,sesid));
        tensor(:,:,trialidx)    = tensor(:,:,trialidx) - repmat(mean(tensor(:,:,trialidx),3),1,1,sum(trialidx));
    end
    
    %Compute noise corr:
    for iCond = 1:nConditions
        trialidx                = splits{iCond}(strcmp(trialData.session_ID,sesid));
        if sum(trialidx)>params.minNtrialscond
            
            for iNx = 1:nNeurons
                for iNy = 1:nNeurons
                    for iT = 1:params.nTimebins
                        noisecorrmat(iCond,neuroncounter + iNx,neuroncounter + iNy,iT) = corr(squeeze(tensor(iNx,iT,trialidx)),squeeze(tensor(iNy,iT,trialidx)));
                    end
                end
            end
            
        end
    end
    
    neuroncounter = neuroncounter + length(tempspikeData.session_ID);
end

%% Set all autocorrelation values of 1 to NaN
noisecorrmat(noisecorrmat>0.99)   = NaN;

%% Report total average noise correlation values during baseline: 

timeidx = params.xtime>=-0.5e6 & params.xtime<=0e6; %Baseline time index

fprintf('Average noise correlation during baseline:\n')
for iCond = 1:nConditions
    temp                        = noisecorrmat(iCond,:,:,timeidx);
    fprintf('%s: %1.4f\n',params.labels_splits{iCond},nanmean(temp(:)))
end

temp                        = nanmean(noisecorrmat(:,:,:,timeidx),4); %Compute noise corr during baseline time window
fprintf('Total average noise correlation during baseline: %1.4f, +- %1.4f std\n',nanmean(temp(:)),nanstd(temp(:))) %report estimate and variability

%% For the statistics (multilevel model):
uMice           = unique(sessionData.mousename);
nMice           = length(uMice);
nNeurons        = length(spikeData.session_ID);
G_mou           = cell(nNeurons,1);
for iMouse = 1:nMice
    G_mou(ismember(spikeData.session_ID,sessionData.session_ID(strcmp(sessionData.mousename,uMice{iMouse})))) = uMice(iMouse);
end

%% Figure 7c: change in noise correlations for visual trial:
timeidx1 = params.xtime>-0.5e6 & params.xtime<=0e6;
timeidx2 = params.xtime>0.2e6  & params.xtime<=1e6;

condidx = [3 4]; %Only visual hits and misses:

figure; set(gcf,'color','w','units','normalized','Position',[0.05 0.05 0.16 0.4])
title('Noise correlations')
hold all;

counter = 0;
           
for iExp = 1:3 %For each experiment:
    idx_exp             = ismember(spikeData.session_ID,sessionData.session_ID(ismember(sessionData.Experiment,params.Experiments(iExp)))); %get all neurons in this task context

    for iCond = condidx %For each quantile reshape spike count correlation values:
        counter                     = counter + 1;
        data_bsl                    = squeeze(nanmean(noisecorrmat(iCond,idx_exp,idx_exp,timeidx1),4));
        data_resp                   = squeeze(nanmean(noisecorrmat(iCond,idx_exp,idx_exp,timeidx2),4));

        %subtract baseline:
        meantoplot(counter,1)    = nanmean(data_resp(:) - data_bsl(:));
        errortoplot(counter,1)   = nanstd(data_resp(:) - data_bsl(:)) / sqrt(sum(~isnan(data_resp(:))));
        
        h = bar(counter,meantoplot(counter,1),0.9);
        h.FaceColor = params.colors_experiments{iExp}; h.EdgeColor = [1 1 1];
        errorbar(counter,meantoplot(counter,1),errortoplot(counter,1),'k.','LineWidth',2)
        
        Y_bsl = [];
        Y_rsp = [];
        G_mou_2 = cell(0,0);
        
        for iMou = 1:nMice
            idx = strcmp(G_mou(idx_exp),uMice(iMou));
            temp = data_bsl(idx,idx);
            Y_bsl = [Y_bsl; temp(:)]; %#ok<*AGROW>
            temp = data_resp(idx,idx);
            Y_rsp = [Y_rsp; temp(:)];
            G_mou_2 = [G_mou_2; repmat(uMice(iMou),numel(temp),1)];
        end
        
        if numel(Y_bsl)>0 && ~all(isnan(Y_bsl))
            Y               = [Y_bsl; Y_rsp];
            X               = [ones(size(Y_bsl)); ones(size(Y_rsp))*2];
            tbl             = table(Y,X,[G_mou_2; G_mou_2],'VariableNames',{'NC','Window','Mouse'}); %Create table for mixed model
            
            lme             = fitlme(tbl,'NC~Window+(1|Mouse)'); %construct linear mixed effects model with fixed effect of temporal window and random intercept for different mice
            stats           = dataset2table(anova(lme,'DFMethod','Satterthwaite')); %Perform ANOVA on model and output as matrix
            fprintf('%s, %s, n=%5.0f pairs: \n',params.ExperimentLabels{iExp}, params.labels_splits{iCond}, sum(~isnan(Y_bsl)))
            fprintf('F(%d,%2.0f)=%1.2f, p=%1.3e; \n',stats{2,3},stats{2,4},stats{2,2},stats{2,5})
            p = stats{2,5};
            if p<0.05
                sigstar([counter-0.15 counter+0.15],p)
            end
            
            tempfile = fullfile(sprintf('SourceData_Fig7c_NC_HitMiss_%s_%s.xlsx',params.ExperimentLabels{iExp}, params.labels_splits{iCond}));
            writetable(tbl,tempfile)

        end
    end
    
end

set(gca,'XTick',1:6,'XTickLabel',repmat(params.labels_splits(condidx),1,3),'XTickLabelRotation',45,'Fontsize',20);
xlim([0.5 2*params.nExperiments+0.5])
ylabel('Delta Rsc','FontSize',20)
ylim([-0.04 0.02])
set(gca,'YTick',[-0.04 -0.02 0 0.02]);


%% Figure S9a: change in noise correlations for visual trial: (full distribution)
figure; set(gcf,'color','w','units','normalized','Position',[0.05 0.05 0.16 0.4])
title('Noise correlations')
hold all;

counter = 0;

for iExp = 1:3 %For each experiment:
    idx_exp             = ismember(spikeData.session_ID,sessionData.session_ID(ismember(sessionData.Experiment,params.Experiments(iExp)))); %get all neurons in this task context
    
    for iCond = condidx %For each quantile reshape spike count correlation values:
        counter                     = counter + 1;
        data_bsl                    = squeeze(nanmean(noisecorrmat(iCond,idx_exp,idx_exp,timeidx1),4));
        data_resp                   = squeeze(nanmean(noisecorrmat(iCond,idx_exp,idx_exp,timeidx2),4));
        
        datatoplot = data_resp(:) - data_bsl(:);
        h = boxplot(datatoplot,'positions',counter,'outliersize',0.0001,'widths',0.5,'whisker',0.5);
        
        set(h(1,1),'Color',params.colors_experiments{iExp},'LineStyle','-','LineWidth',1);
        set(h(2,1),'Color',params.colors_experiments{iExp},'LineStyle','-','LineWidth',1);
        set(h(3,1),'Color',params.colors_experiments{iExp});
        set(h(4,1),'Color',params.colors_experiments{iExp});
        set(h(5,1),'Color',params.colors_experiments{iExp},'LineWidth',1);
        set(h(6,1),'Color',params.colors_experiments{iExp},'LineWidth',1);
        delete(h(7,1));
    end
end

set(gca,'XTick',1:6,'XTickLabel',repmat(params.labels_splits(condidx),1,3),'XTickLabelRotation',45,'Fontsize',20);
plot([-1 10],[0 0],'k-','LineWidth',0.5)
xlim([0.5 2*params.nExperiments+0.5])
ylabel('Delta Rsc','FontSize',20)
ylim([-0.275 0.275])
set(gca,'YTick',-0.5:0.1:0.5);

%% Figure S9b: change in noise correlations for auditory trials:
timeidx1 = params.xtime>-0.5e6 & params.xtime<=0e6;
timeidx2 = params.xtime>0.2e6  & params.xtime<=1e6;

condidx = [1 2]; %Only auditory hits and misses:

figure; set(gcf,'color','w','units','normalized','Position',[0.05 0.05 0.16 0.4])
title('Noise correlations')
hold all;

counter = 0;
           
for iExp = [1 3] %For each experiment:
    idx_exp             = ismember(spikeData.session_ID,sessionData.session_ID(ismember(sessionData.Experiment,params.Experiments(iExp)))); %get all neurons in this task context

    for iCond = condidx %For each quantile reshape spike count correlation values:
        counter                     = counter + 1;
        data_bsl                    = squeeze(nanmean(noisecorrmat(iCond,idx_exp,idx_exp,timeidx1),4));
        data_resp                   = squeeze(nanmean(noisecorrmat(iCond,idx_exp,idx_exp,timeidx2),4));

        %subtract baseline:
        meantoplot(counter,1)    = nanmean(data_resp(:) - data_bsl(:));
        errortoplot(counter,1)   = nanstd(data_resp(:) - data_bsl(:)) / sqrt(sum(~isnan(data_resp(:))));
        
        h = bar(counter,meantoplot(counter,1),0.9);
        h.FaceColor = params.colors_experiments{iExp}; h.EdgeColor = [1 1 1];
        errorbar(counter,meantoplot(counter,1),errortoplot(counter,1),'k.','LineWidth',2)
        
        Y_bsl = [];
        Y_rsp = [];
        G_mou_2 = cell(0,0);
        
        for iMou = 1:nMice
            idx = strcmp(G_mou(idx_exp),uMice(iMou));
            temp = data_bsl(idx,idx);
            Y_bsl = [Y_bsl; temp(:)]; %#ok<*AGROW>
            temp = data_resp(idx,idx);
            Y_rsp = [Y_rsp; temp(:)];
            G_mou_2 = [G_mou_2; repmat(uMice(iMou),numel(temp),1)];
        end
        
        if numel(Y_bsl)>0 && ~all(isnan(Y_bsl))
            Y               = [Y_bsl; Y_rsp];
            X               = [ones(size(Y_bsl)); ones(size(Y_rsp))*2];
            tbl             = table(Y,X,[G_mou_2; G_mou_2],'VariableNames',{'NC','Window','Mouse'}); %Create table for mixed model
            
            lme             = fitlme(tbl,'NC~Window+(1|Mouse)'); %construct linear mixed effects model with fixed effect of temporal window and random intercept for different mice
            stats           = dataset2table(anova(lme,'DFMethod','Satterthwaite')); %Perform ANOVA on model and output as matrix
            fprintf('%s, %s, n=%5.0f pairs: \n',params.ExperimentLabels{iExp}, params.labels_splits{iCond}, sum(~isnan(Y_bsl)))
            fprintf('F(%d,%2.0f)=%1.2f, p=%1.3e; \n',stats{2,3},stats{2,4},stats{2,2},stats{2,5})
            p = stats{2,5};
            if p<0.05
                sigstar([counter-0.15 counter+0.15],p)
            end
        end
    end
end

set(gca,'XTick',1:6,'XTickLabel',repmat(params.labels_splits(condidx),1,3),'XTickLabelRotation',45,'Fontsize',20);
xlim([0.5 2*params.nExperiments+0.5])
ylabel('Delta Rsc','FontSize',20)
ylim([-0.04 0.02])
set(gca,'YTick',[-0.04 -0.02 0 0.02]);


%% Figure 9b: change in noise correlations for visual trial: (full distribution)
figure; set(gcf,'color','w','units','normalized','Position',[0.05 0.05 0.16 0.4])
title('Noise correlations')
hold all;

counter = 0;

for iExp = [1 3] %For each experiment:
    idx_exp             = ismember(spikeData.session_ID,sessionData.session_ID(ismember(sessionData.Experiment,params.Experiments(iExp)))); %get all neurons in this task context
    
    for iCond = condidx %For each quantile reshape spike count correlation values:
        counter                     = counter + 1;
        data_bsl                    = squeeze(nanmean(noisecorrmat(iCond,idx_exp,idx_exp,timeidx1),4));
        data_resp                   = squeeze(nanmean(noisecorrmat(iCond,idx_exp,idx_exp,timeidx2),4));
        
        datatoplot = data_resp(:) - data_bsl(:);
        h = boxplot(datatoplot,'positions',counter,'outliersize',0.0001,'widths',0.5,'whisker',0.5);
        
        set(h(1,1),'Color',params.colors_experiments{iExp},'LineStyle','-','LineWidth',1);
        set(h(2,1),'Color',params.colors_experiments{iExp},'LineStyle','-','LineWidth',1);
        set(h(3,1),'Color',params.colors_experiments{iExp});
        set(h(4,1),'Color',params.colors_experiments{iExp});
        set(h(5,1),'Color',params.colors_experiments{iExp},'LineWidth',1);
        set(h(6,1),'Color',params.colors_experiments{iExp},'LineWidth',1);
        delete(h(7,1));
    end
end

set(gca,'XTick',1:6,'XTickLabel',repmat(params.labels_splits(condidx),1,3),'XTickLabelRotation',45,'Fontsize',20);
xlim([0.5 2*params.nExperiments+0.5])
ylabel('Delta Rsc','FontSize',20)
ylim([-0.2 0.2])
plot([-1 10],[0 0],'k-','LineWidth',0.5)
set(gca,'YTick',-0.5:0.1:0.5);

%% Figure: noise correlations over time (aligned to stimulus): (not in MS)
params.baselinecorr     = 1;
params.idx_bsl = params.xtime<0;

for iExp = 1:3 %For each task contingency cohort:
    idx_exp             = ismember(spikeData.session_ID,sessionData.session_ID(ismember(sessionData.Experiment,params.Experiments(iExp)))); %get all neurons in this task context
    
    nCross              = sum(idx_exp)^2; %calc max combinations
    datatoplot          = NaN(nConditions,nCross,params.nTimebins); %init var
    
    for iCond = 1:nConditions %For each condition reshape spike count correlation values:
        temp                        = noisecorrmat(iCond,idx_exp,idx_exp,:);
        temp                        = reshape(temp,sum(idx_exp)^2,params.nTimebins);
        datatoplot(iCond,:,:)       = temp;
    end
    
    idx             = any(~isnan(datatoplot),3); %calculate how many correlation values:
    fprintf('n = %d pairwise combinations\n',sum(sum(idx)))
    
    %Compute mean and sem:
    meantoplot      = squeeze(nanmean(datatoplot,2));
    errortoplot     = squeeze(nanstd(datatoplot,[],2) ./ sqrt(sum(~isnan(datatoplot),2)));
    
    if params.baselinecorr %Subtract baseline values 
        meantoplot      = meantoplot - repmat(nanmean(meantoplot(:,params.idx_bsl),2),1,params.nTimebins);
    end
    
    %Figure of overall averaged noise correlations for the three reaction time quantiles:
    figure; hold all; set(gcf,'units','normalized','Position',[0.45 0.45 0.32 0.3],'color','w');
    handles = [];
    
    for iCond = 1:nConditions
        h = shadedErrorBar(params.xtime, meantoplot(iCond,:),errortoplot(iCond,:),{params.lines_splits{iCond},'LineWidth',3,'Color',params.colors_trialtypes{ceil(iCond/2)}},1);
        handles(end+1) = h.mainLine; %#ok<SAGROW>
        delete(h.edge(1)); delete(h.edge(2));
        %         handles(end+1) = plot(params.xtime, meantoplot(iCond,:) - nanmean(meantoplot(iCond,params.xtime<0)),params.lines_splits{iCond},'LineWidth',4,'Color',params.colors_splits(iCond,:));
    end
        
    %Figure makeup:
    plot([0 0],[-0.2 0.2],'k','LineWidth',2)
    ylabel('Baseline corrected rsc','FontSize',20)
    legend(handles,params.labels_splits); legend boxoff;
    xlim([-0.2e6 1.2e6])
    ylim([0 0.13])
    if params.baselinecorr %Subtract baseline values from
        ylim([-0.07 0.05])
    end
    xlabel('Time (ms)')
    set(gca,'XTick',-2e6:0.2e6:2e6,'XTickLabels',(-2e6:0.2e6:2e6)*1e-3,'YTick',[-0.1 -0.05 0 0.05 0.1])
    grid on;
    
    %     export_fig(gcf,fullfile(params.savedir,sprintf('V1_NC_trialtypes_%s_alignStim.eps',params.ExperimentLabels{iExp})))
end

%%















%% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Figure 2:
% Show noise correlations over time for each of quantiles of reaction times
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%% Labels and colors:
params.nQuantiles                  = 3;
params.colors_quantiles            = hot(5*params.nQuantiles);
params.colors_quantiles            = params.colors_quantiles(2:3:2+3*(params.nQuantiles-1),:);
params.colors_quantiles            = [0.5 0.5 0.5; params.colors_quantiles;];
params.lines_quantiles             = {'-.' '-' '-' '-'};
params.labels_quantiles            = {'Catch' 'Fast' 'Mid' 'Slow'};

%% Trial splits based on reaction time:

%Parameters for slow/fast reaction time splits:
RT_quantile = NaN(2,3);
figure;  hold all; set(gcf,'units','normalized','Position',[0.45 0.05 0.3 0.4],'color','w');

for iExp = [2 3]
    %Determine RT quantiles for each experiment separately:
    [~,temptrialData]                       = MOL_getTempPerSes(sessionData.session_ID(ismember(sessionData.Experiment,params.Experiments(iExp))),sessionData,trialData);%Get the sessionData for each session individually:
    all_RT                                  = temptrialData.responseLatency(strcmp(temptrialData.trialType,'X') & temptrialData.vecResponse==2);
    
    %Define quantile edges for trials:
    [edges(iExp-1,:)]                       = quantile(all_RT,linspace(0,1,params.nQuantiles+1));
    
    % Show distribution of reaction times
    binres      = 10e3;
    
    binedges    = binres/2:binres:1e6;
    histx       = binedges(1:end-1)+binres/2;
    histy       = histcounts(all_RT,binedges,'normalization','count');
    
    for iQ = 1:params.nQuantiles
        idx = histx>=edges(iExp-1,iQ) & histx<=edges(iExp-1,iQ+1);
        bar(histx(idx),-(iExp*2-5)*histy(idx),'FaceColor',params.colors_quantiles(iQ+1,:),'EdgeColor','k')
        
        RT_quantile(iExp-1,iQ) = nanmedian(all_RT(all_RT>=edges(iExp-1,iQ) & all_RT<=edges(iExp-1,iQ+1)));
    end
        
end

text(700e3,60,'UST','FontSize',15)
text(700e3,-60,'MST','FontSize',15)

xlim([200e3 1e6])
ylim([-80 80])
set(gca,'XTick',(200:200:1000)*1e3,'XTickLabels',200:200:1000,'YTick',[-80 80])
xlabel('Reaction time (ms)')
ylabel('Count')
legend(params.labels_quantiles(2:end)); legend boxoff;


%% Statistics: comparing RT for UST and MST:
% iExp = 2;
% [~,temptrialData]                           = MOL_getTempPerSes(sessionData.session_ID(ismember(sessionData.Experiment,params.Experiments(iExp))),sessionData,trialData);%Get the sessionData for each session individually:
% all_RT_UST                                  = temptrialData.responseLatency(strcmp(temptrialData.trialType,'X') & temptrialData.vecResponse==2);
% iExp = 3;
% [~,temptrialData]                           = MOL_getTempPerSes(sessionData.session_ID(ismember(sessionData.Experiment,params.Experiments(iExp))),sessionData,trialData);%Get the sessionData for each session individually:
% all_RT_MST                                  = temptrialData.responseLatency(strcmp(temptrialData.trialType,'X') & temptrialData.vecResponse==2);
% 
% p = ranksum(all_RT_UST,all_RT_MST);
% fprintf('%d UST visual hits vs %d MST visual hits: p=%1.3e\n',numel(all_RT_UST),numel(all_RT_MST),p)

%% Statistics: comparing RT for UST and MST:
nSessions                   = length(sessionData.session_ID);
all_RT = NaN(nSessions,1);
exp = cell(nSessions,1);

for iSes = 1:nSessions
    [~,temptrialData]                           = MOL_getTempPerSes(sessionData.session_ID(iSes),sessionData,trialData);%Get the sessionData for each session individually:
    all_RT(iSes)                                = nanmedian(temptrialData.responseLatency(strcmp(temptrialData.trialType,'X') & temptrialData.vecResponse==2));
    exp(iSes)                                   = sessionData.Experiment(iSes);
end

p = ranksum(all_RT(strcmp(exp,params.Experiments{2})),all_RT(strcmp(exp,params.Experiments{3})));
fprintf('%d UST sessions vs %d MST sessions: p=%1.4f\n',sum(strcmp(exp,params.Experiments{2})),sum(strcmp(exp,params.Experiments{3})),p)

tbl             = table(all_RT,exp,'VariableNames',{'Median RT (visual trials)','Cohort'}); %Create table for mixed model
tempfile        = fullfile('SourceData_Fig7d_RTsessions_Cohort.xlsx');
writetable(tbl,tempfile)
    
%% Get indices of the trials per quantile:
splits                      = {};
splits{end+1}               = strcmp(trialData.trialType,'P') & trialData.vecResponse==3;
splits{end+1}               = false(size(trialData.trialType));
splits{end+1}               = false(size(trialData.trialType));
splits{end+1}               = false(size(trialData.trialType));

for iExp = [2 3]
    idx_exp                 = ismember(trialData.session_ID,sessionData.session_ID(ismember(sessionData.Experiment,params.Experiments(iExp))));
    for iQ = 1:params.nQuantiles
        idx_qua             = strcmp(trialData.trialType,'X') & trialData.vecResponse==2 & (trialData.responseLatency>=edges(iExp-1,iQ) & trialData.responseLatency<=edges(iExp-1,iQ+1));
        splits{iQ+1}        = splits{iQ+1} | (idx_qua & idx_exp);
    end
end

%% Compute noise correlations:

%Some measures of dimensionality:
nTotalNeurons               = length(spikeData.session_ID);
nSessions                   = length(sessionData.session_ID);

params                      = params_histresponse_coding(params);
params.binsize              = 25e3;
params.smoothing            = 1;
params.conv_win             = 'chg';       %Type of window used for smoothing {flat, gaussian)
params.conv_sigma           = 0.100e6;           %sd of gaussian window for smoothing
params.conv_twin            = round((9*params.conv_sigma)/params.binsize)*params.binsize;         %Window size for smoothing

%Construct bin edges and time axis
params.edges                = [params.t_pre:params.binsize:params.t_post] - params.binsize/2;                    %#ok<NBRAK> Define edges of the bins
params.xtime                = params.t_pre:params.binsize:params.t_post-params.binsize;                        %Define time axis
params.nTimebins            = length(params.xtime); %number of time bins

params.AlignOn              = 'stimChange';

params.minNtrialscond       = 10;

%Construct bin edges and time axis
params.edges                = [params.t_pre:params.binsize:params.t_post] - params.binsize/2;                    %#ok<NBRAK> Define edges of the bins
params.xtime                = params.t_pre:params.binsize:params.t_post-params.binsize;                        %Define time axis
params.nTimebins            = length(params.xtime); %number of time bins

nConditions                 = numel(splits);
noisecorrmat                = NaN(nConditions,nTotalNeurons,nTotalNeurons,params.nTimebins); %init output matrix

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
    
    %subtract mean response for each condition:
    for iCond = 1:nConditions
        trialidx                = splits{iCond}(strcmp(trialData.session_ID,sesid));
        tensor(:,:,trialidx)    = tensor(:,:,trialidx) - repmat(nanmean(tensor(:,:,trialidx),3),1,1,sum(trialidx));
    end
    
    %Compute noise corr:
    for iCond = 1:nConditions
        trialidx                = splits{iCond}(strcmp(trialData.session_ID,sesid));
        if sum(trialidx)>params.minNtrialscond
            
            for iNx = 1:nNeurons
                for iNy = 1:nNeurons
                    for iT = 1:params.nTimebins
                        noisecorrmat(iCond,neuroncounter + iNx,neuroncounter + iNy,iT) = corr(squeeze(tensor(iNx,iT,trialidx)),squeeze(tensor(iNy,iT,trialidx)));
                    end
                end
            end
            
        end
    end
    
    neuroncounter = neuroncounter + length(tempspikeData.session_ID);
end

%% Set all autocorrelation values of 1 to NaN
noisecorrmat(noisecorrmat>0.99)   = NaN;

%%
params.baselinecorr     = 1;
params.nstd             = 2;
params.idx_bsl          = params.xtime>=-.5e6 & params.xtime<=0e6;

tBin                    = NaN(2,nConditions);

for iExp = [2 3] %Only for trained task contexts:
    idx_exp             = ismember(spikeData.session_ID,sessionData.session_ID(ismember(sessionData.Experiment,params.Experiments(iExp)))); %get all neurons in this task context
    idx_exp             = idx_exp & spikeData.avg_fr>1; %get all neurons in this task context
    nCross              = sum(idx_exp)^2; %calc max combinations
    datatoplot          = NaN(nConditions,nCross,params.nTimebins); %init var
    
    for iCond = 1:nConditions %For each quantile reshape spike count correlation values:
        temp                        = noisecorrmat(iCond,idx_exp,idx_exp,:);
        temp                        = reshape(temp,sum(idx_exp)^2,params.nTimebins);
        datatoplot(iCond,:,:)       = temp;
    end
    
    idx             = any(~isnan(datatoplot),3); %calculate how many correlation values:
    fprintf('n = %d pairwise combinations\n',sum(sum(idx)))
    
    %Compute mean and sem:
    meantoplot      = squeeze(nanmean(datatoplot,2));
    errortoplot     = squeeze(nanstd(datatoplot,[],2) ./ sqrt(sum(~isnan(datatoplot),2)));
    
    if params.baselinecorr %Subtract baseline values from 
        meantoplot      = meantoplot - repmat(nanmean(meantoplot(:,params.idx_bsl),2),1,params.nTimebins);
    end
    
    %Figure of overall averaged noise correlations for the three reaction time quantiles:
    figure; hold all; set(gcf,'units','normalized','Position',[0.45 0.45 0.32 0.3],'color','w');
    handles = [];
    
    for iCond = 2:nConditions
        h = shadedErrorBar(params.xtime, meantoplot(iCond,:),errortoplot(iCond,:),{params.lines_quantiles{iCond},'LineWidth',3,'Color',params.colors_quantiles(iCond,:)},0);
        handles(end+1) = h.mainLine; %#ok<SAGROW>
        delete(h.edge(1)); delete(h.edge(2)); %h.edge(1).Color = [1 1 1];  h.edge(2).Color =  [1 1 1];
    end
        
    for iCond = 2:nConditions %Skip catch, start at 2
        
        temp = NaN(1000,1);
        tempidx = find(~isnan(datatoplot(iCond,:,1)));
        if tempidx
           
            CImean = repmat(nanmean(meantoplot(iCond,params.idx_bsl)),1,params.nTimebins);
            CIband = repmat(params.nstd*nanstd(meantoplot(iCond,params.idx_bsl)),1,params.nTimebins);
           
            plot(params.xtime,CImean-CIband,'--','LineWidth',1,'Color',params.colors_quantiles(iCond,:));
            
            temp = find(meantoplot(iCond,:)<(CImean-CIband) & params.xtime>0,1);
            if temp
                tBin(iExp-1,iCond) = temp;
                plot([params.xtime(tBin(iExp-1,iCond)) params.xtime(tBin(iExp-1,iCond))],[meantoplot(iCond,tBin(iExp-1,iCond)) 0.08],'Color',params.colors_quantiles(iCond,:))
                %             plot([params.xtime(tBin(iCond)) params.xtime(tBin(iCond))],[meantoplot(iCond,tBin(iCond)) 0.08],'Color',params.colors_trialtypes{ceil(iCond/4)})
            end
        end
    end
    for iCond = 1:nConditions
        if ~isnan(tBin(iExp-1,iCond))
            fprintf('Earliest timepoint of significant decorrelation for %s, %s: %3.0f ms\n',params.ExperimentLabels{iExp},params.labels_quantiles{iCond},params.xtime(tBin(iExp-1,iCond))*1e-3)
        end
    end
    
    %Figure makeup:
    plot([0 0],[-0.2 0.2],'k','LineWidth',2)
    ylabel('Baseline corrected rsc','FontSize',20)
    legend(handles,params.labels_quantiles(2:end)); legend boxoff;
    xlim([-0.3e6 1e6])
    if params.baselinecorr
        ylim([-0.075 0.05])
    else
        ylim([-0.02 0.12])
    end
    xlabel('Time (ms)')
    set(gca,'XTick',-1.2e6:0.2e6:1e6,'XTickLabels',(-1.2e6:0.2e6:1e6)*1e-3,'YTick',[-0.1 -0.075 0 0.05 0.1])
%     grid on;
    
    %     export_fig(gcf,fullfile(params.savedir,sprintf('V1_NC_Quantiles_%s_alignStim.eps',params.ExperimentLabels{iExp})))
end

%% Correlate response latency with onset of decorrelation:

RT_toplot       = RT_quantile;
tDecor_toplot   = params.xtime(tBin(:,2:end));

figure; hold all; set(gcf,'units','normalized','Position',[0.45 0.45 0.16 0.23],'color','w');
scatter(RT_toplot(1,:)*1e-3,tDecor_toplot(1,:)*1e-3,120,params.colors_quantiles(2:end,:),'o','filled')
scatter(RT_toplot(2,:)*1e-3,tDecor_toplot(2,:)*1e-3,120,params.colors_quantiles(2:end,:),'s','filled')
plot([0 1e6],[0 1e6],'k:','LineWidth',2)
xlabel('Median RT')
ylabel('Earliest Decor')
set(gca,'XTick',[0 500 1000],'YTick',[0 500 1000])
xlim([0 1000])
ylim([0 1000])

[r,p] = corr(RT_toplot(:),tDecor_toplot(:));
text(100,600,sprintf('r=%1.3f,p=%1.3f',r,p),'FontSize',12)
legend({'UST' 'MST'},'Location','NorthWest'); legend boxoff

tempfile = fullfile('SourceData_Fig7f_RT_tDecor.xlsx');
writetable(table(RT_toplot(:),tDecor_toplot(:),'VariableNames',{'Reaction Time','t_Decorrelation'}),tempfile)

%%









%% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Figure 3:
% Show noise correlations over time for each of quantiles of reaction times but ALIGNED TO FIRST LICK
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%% Compute noise correlations:

%Create timestamp of the first lick:
trialData.firstLick         = trialData.stimChange;
trialData.firstLick(~isnan(trialData.responseLatency)) = trialData.firstLick(~isnan(trialData.responseLatency)) + trialData.responseLatency(~isnan(trialData.responseLatency));

% add pseudo reaction times for no lick trials:
z = trialData.responseLatency(trialData.responseLatency>150e3);
h = fitdist(z,'lognormal');

nNoLick                     = sum(isnan(trialData.responseLatency));
%lognrnd generates a random number from the lognormal distribution with the distribution parameters mu (mean of logarithmic values) and sigma (standard deviation of logarithmic values).
newRT                       = lognrnd(h.mu,h.sigma,[nNoLick 1]);
newRT(newRT<100e3)          = lognrnd(h.mu,h.sigma,[sum(newRT<100e3) 1]);
trialData.firstLick(isnan(trialData.responseLatency)) = trialData.stimChange(isnan(trialData.responseLatency)) + newRT;


%% To show fit: Uncomment to display:
% figure; hold all;
% edges = 100e3:20e3:1.5e6;
% x = histcounts(trialData.responseLatency,edges,'Normalization','probability');
% plot(edges(2:end),x/max(x),'k.','MarkerSize',25);
% fitpdf = pdf(h,edges);
% plot(edges,fitpdf/max(fitpdf),'r-','LineWidth',2);

%% Compute noise correlations:
params.AlignOn              = 'firstLick';

params.t_pre                = -2.5e6;
params.t_post               = +2.5e6;

%Construct bin edges and time axis
params.edges                = [params.t_pre:params.binsize:params.t_post] - params.binsize/2;                    %#ok<NBRAK> Define edges of the bins
params.xtime                = params.t_pre:params.binsize:params.t_post-params.binsize;                        %Define time axis
params.nTimebins            = length(params.xtime); %number of time bins


noisecorrmat                = NaN(nConditions,nTotalNeurons,nTotalNeurons,params.nTimebins); %init output matrix

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
    
    %subtract mean response for each condition:
    for iCond = 1:nConditions
        trialidx                = splits{iCond}(strcmp(trialData.session_ID,sesid));
        tensor(:,:,trialidx)    = tensor(:,:,trialidx) - repmat(mean(tensor(:,:,trialidx),3),1,1,sum(trialidx));
    end
    
    %Compute noise corr:
    for iCond = 1:nConditions
        trialidx                = splits{iCond}(strcmp(trialData.session_ID,sesid));
        if sum(trialidx)>params.minNtrialscond
            
            for iNx = 1:nNeurons
                for iNy = 1:nNeurons
                    for iT = 1:params.nTimebins
                        noisecorrmat(iCond,neuroncounter + iNx,neuroncounter + iNy,iT) = corr(squeeze(tensor(iNx,iT,trialidx)),squeeze(tensor(iNy,iT,trialidx)));
                    end
                end
            end
            
        end
    end
    
    neuroncounter = neuroncounter + length(tempspikeData.session_ID);
end

%% Set all autocorrelation values of 1 to NaN
noisecorrmat(noisecorrmat>0.99)   = NaN;

%% Make figure:
params.baselinecorr     = 1;
params.idx_bsl          = params.xtime>-1e6 & params.xtime<-0.3e6;
params.bsl_end          = params.xtime(find(params.idx_bsl,1,'last'));

params.nstd             = 2;

tBin                    = NaN(2,nConditions);

for iExp = [2 3] %Only for trained task contexts:
    idx_exp             = ismember(spikeData.session_ID,sessionData.session_ID(ismember(sessionData.Experiment,params.Experiments(iExp)))); %get all neurons in this task context
    idx_exp             = idx_exp & spikeData.avg_fr>1; %get all neurons in this task context
    nCross              = sum(idx_exp)^2; %calc max combinations
    datatoplot          = NaN(nConditions,nCross,params.nTimebins); %init var
    
    for iCond = 1:nConditions %For each quantile reshape spike count correlation values:
        temp                        = noisecorrmat(iCond,idx_exp,idx_exp,:);
        temp                        = reshape(temp,sum(idx_exp)^2,params.nTimebins);
        datatoplot(iCond,:,:)       = temp;
    end
    
    idx             = any(~isnan(datatoplot),3); %calculate how many correlation values:
    fprintf('n = %d pairwise combinations\n',sum(sum(idx)))
    
    %Compute mean and sem:
    meantoplot      = squeeze(nanmean(datatoplot,2));
    errortoplot     = squeeze(nanstd(datatoplot,[],2) ./ sqrt(sum(~isnan(datatoplot),2)));
    
    if params.baselinecorr %Subtract baseline values from 
        meantoplot      = meantoplot - repmat(nanmean(meantoplot(:,params.idx_bsl),2),1,params.nTimebins);
    end
    
    %Figure of overall averaged noise correlations for the three reaction time quantiles:
    figure; hold all; set(gcf,'units','normalized','Position',[0.45 0.45 0.32 0.3],'color','w');
    handles = [];
    
    for iCond = 2:nConditions
        h = shadedErrorBar(params.xtime, meantoplot(iCond,:),errortoplot(iCond,:),{params.lines_quantiles{iCond},'LineWidth',3,'Color',params.colors_quantiles(iCond,:)},0);
        handles(end+1) = h.mainLine; %#ok<SAGROW>
        delete(h.edge(1)); delete(h.edge(2)); %h.edge(1).Color = [1 1 1];  h.edge(2).Color =  [1 1 1];
    end
        
    for iCond = 2:nConditions
        
        temp        = NaN(1000,1);
        tempidx     = find(~isnan(datatoplot(iCond,:,1)));
        if tempidx
            
            CImean = repmat(nanmean(meantoplot(iCond,params.idx_bsl)),1,params.nTimebins);
            CIband = repmat(params.nstd*nanstd(meantoplot(iCond,params.idx_bsl)),1,params.nTimebins);
            
            h = plot(params.xtime,CImean-CIband,'--','LineWidth',1,'Color',params.colors_quantiles(iCond,:));
            
            temp = find(meantoplot(iCond,:)<(CImean-CIband) & params.xtime>params.bsl_end,1);
            if temp
                tBin(iExp-1,iCond) = temp;
                plot([params.xtime(tBin(iExp-1,iCond)) params.xtime(tBin(iExp-1,iCond))],[meantoplot(iCond,tBin(iExp-1,iCond)) 0.08],'Color',params.colors_quantiles(iCond,:))
            end
        end
    end
    for iCond = 1:nConditions
        if ~isnan(tBin(iExp-1,iCond))
            fprintf('Earliest timepoint of significant decorrelation for %s, %s: %3.0f ms\n',params.ExperimentLabels{iExp},params.labels_quantiles{iCond},params.xtime(tBin(iExp-1,iCond))*1e-3)
        end
    end
    
    %Figure makeup:
    plot([0 0],[-0.2 0.2],'k','LineWidth',2)
    ylabel('Baseline corrected rsc','FontSize',20)
    legend(handles,params.labels_quantiles(2:end)); legend boxoff;
    xlim([-0.6e6 1e6])
    ylim([-0.09 0.05])
    xlabel('Time (ms)')
    set(gca,'XTick',-2e6:0.2e6:2e6,'XTickLabels',(-2e6:0.2e6:2e6)*1e-3,'YTick',[-0.1 -0.05 0 0.05 0.1])
    grid on;
    
    %     export_fig(gcf,fullfile(params.savedir,sprintf('V1_NC_Quantiles_%s_alignStim.eps',params.ExperimentLabels{iExp})))
end

%% Correlate response latency with onset of decorrelation:

RT_toplot       = RT_quantile;
tDecor_toplot   = params.xtime(tBin(:,2:end));

figure; hold all; set(gcf,'units','normalized','Position',[0.45 0.45 0.16 0.23],'color','w');
scatter(RT_toplot(1,:)*1e-3,tDecor_toplot(1,:)*1e-3,120,params.colors_quantiles(2:end,:),'o','filled')
scatter(RT_toplot(2,:)*1e-3,tDecor_toplot(2,:)*1e-3,120,params.colors_quantiles(2:end,:),'s','filled')
plot([0 1e6],[0 1e6],'k:','LineWidth',2)
xlabel('Median RT')
ylabel('Earliest Decor')
set(gca,'XTick',[0 500 1000],'YTick',[-500 -250 0])
xlim([0 1000])
ylim([-500 0])

[r,p] = corr(RT_toplot(:),tDecor_toplot(:));
text(100,-400,sprintf('r=%1.3f,p=%1.3f',r,p),'FontSize',12)
legend({'UST' 'MST'},'Location','NorthWest'); legend boxoff

tempfile = fullfile('SourceData_FigS9c_RT_DecorLickAlign.xlsx');
writetable(table(RT_toplot(:),tDecor_toplot(:),'VariableNames',{'RT' 'tDecorrelation'}),tempfile)

%%















%% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Figure supplementary: CONTROL ANALYSIS WITH CATCH TRIALS:
% Show noise correlations over time for each of quantiles of reaction times for catch trials
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%% Labels and colors:
params.nQuantiles                  = 3;
params.colors_quantiles            = hot(5*params.nQuantiles);
params.colors_quantiles            = params.colors_quantiles(2:3:2+3*(params.nQuantiles-1),:);
params.colors_quantiles            = [0.5 0.5 0.5; params.colors_quantiles;];
params.lines_quantiles             = {'-.' '-' '-' '-'};
params.labels_quantiles            = {'Catch' 'Fast' 'Mid' 'Slow'};

%% Trial splits based on reaction time:

%Parameters for slow/fast reaction time splits:
RT_quantile = NaN(2,3);
figure;  hold all; set(gcf,'units','normalized','Position',[0.45 0.05 0.3 0.4],'color','w');

for iExp = [2 3]
    %Determine RT quantiles for each experiment separately:
    [~,temptrialData]                       = MOL_getTempPerSes(sessionData.session_ID(ismember(sessionData.Experiment,params.Experiments(iExp))),sessionData,trialData);%Get the sessionData for each session individually:
    all_RT                                  = temptrialData.responseLatency(strcmp(temptrialData.trialType,'P') & temptrialData.vecResponse==2 & temptrialData.responseLatency>200e3);
%     all_RT                                  = temptrialData.responseLatency(strcmp(temptrialData.trialType,'Y') & temptrialData.vecResponse==1);
    
    %Define quantile edges for trials:
    [edges(iExp-1,:)]                       = quantile(all_RT,linspace(0,1,params.nQuantiles+1));
    
    % Show distribution of reaction times
    binres      = 10e3;
    
    binedges    = binres/2:binres:1e6;
    histx       = binedges(1:end-1)+binres/2;
    histy       = histcounts(all_RT,binedges,'normalization','count');
    
    for iQ = 1:params.nQuantiles
        idx = histx>=edges(iExp-1,iQ) & histx<=edges(iExp-1,iQ+1);
        bar(histx(idx),-(iExp*2-5)*histy(idx),'FaceColor',params.colors_quantiles(iQ+1,:),'EdgeColor','k')
        
        RT_quantile(iExp-1,iQ) = nanmedian(all_RT(all_RT>=edges(iExp-1,iQ) & all_RT<=edges(iExp-1,iQ+1)));
    end
        
end

text(700e3,60,'UST','FontSize',15)
text(700e3,-60,'MST','FontSize',15)

xlim([0e3 1e6])
ylim([-80 80])
set(gca,'XTick',(200:200:1000)*1e3,'XTickLabels',200:200:1000,'YTick',[-80 80])
xlabel('Reaction time (ms)')
ylabel('Count')
legend(params.labels_quantiles(2:end)); legend boxoff;

iExp = 2;
[~,temptrialData]                       = MOL_getTempPerSes(sessionData.session_ID(ismember(sessionData.Experiment,params.Experiments(iExp))),sessionData,trialData);%Get the sessionData for each session individually:
all_RT_UST                                  = temptrialData.responseLatency(strcmp(temptrialData.trialType,'P') & temptrialData.vecResponse==2);
iExp = 3;
[~,temptrialData]                       = MOL_getTempPerSes(sessionData.session_ID(ismember(sessionData.Experiment,params.Experiments(iExp))),sessionData,trialData);%Get the sessionData for each session individually:
all_RT_MST                                  = temptrialData.responseLatency(strcmp(temptrialData.trialType,'P') & temptrialData.vecResponse==2);

p = ranksum(all_RT_UST,all_RT_MST);
fprintf('%d UST false alarms vs %d MST false alarms: p=%1.3e\n',numel(all_RT_UST),numel(all_RT_MST),p)

%% Get indices of the trials per quantile:
splits                      = {};
splits{end+1}               = strcmp(trialData.trialType,'P') & trialData.vecResponse==3;
splits{end+1}               = false(size(trialData.trialType));
splits{end+1}               = false(size(trialData.trialType));
splits{end+1}               = false(size(trialData.trialType));

for iExp = [2 3]
    idx_exp                 = ismember(trialData.session_ID,sessionData.session_ID(ismember(sessionData.Experiment,params.Experiments(iExp))));
    for iQ = 1:params.nQuantiles
        idx_qua             = strcmp(trialData.trialType,'P') & trialData.vecResponse==2 & (trialData.responseLatency>=edges(iExp-1,iQ) & trialData.responseLatency<=edges(iExp-1,iQ+1));
        splits{iQ+1}        = splits{iQ+1} | (idx_qua & idx_exp);
    end
end

%% Compute noise correlations:

%Some measures of dimensionality:
nTotalNeurons               = length(spikeData.session_ID);
nSessions                   = length(sessionData.session_ID);

params                      = params_histresponse_coding(params);
params.binsize              = 25e3;
params.smoothing            = 1;
params.conv_win             = 'chg';       %Type of window used for smoothing {flat, gaussian)
params.conv_sigma           = 0.100e6;           %sd of gaussian window for smoothing
params.conv_twin            = round((9*params.conv_sigma)/params.binsize)*params.binsize;         %Window size for smoothing

%Construct bin edges and time axis
params.edges                = [params.t_pre:params.binsize:params.t_post] - params.binsize/2;                    %#ok<NBRAK> Define edges of the bins
params.xtime                = params.t_pre:params.binsize:params.t_post-params.binsize;                        %Define time axis
params.nTimebins            = length(params.xtime); %number of time bins

params.AlignOn              = 'stimChange';

params.minNtrialscond       = 10;

%Construct bin edges and time axis
params.edges                = [params.t_pre:params.binsize:params.t_post] - params.binsize/2;                    %#ok<NBRAK> Define edges of the bins
params.xtime                = params.t_pre:params.binsize:params.t_post-params.binsize;                        %Define time axis
params.nTimebins            = length(params.xtime); %number of time bins

nConditions                 = numel(splits);
noisecorrmat                = NaN(nConditions,nTotalNeurons,nTotalNeurons,params.nTimebins); %init output matrix

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
    
    %subtract mean response for each condition:
    for iCond = 1:nConditions
        trialidx                = splits{iCond}(strcmp(trialData.session_ID,sesid));
        tensor(:,:,trialidx)    = tensor(:,:,trialidx) - repmat(nanmean(tensor(:,:,trialidx),3),1,1,sum(trialidx));
    end
    
    %Compute noise corr:
    for iCond = 1:nConditions
        trialidx                = splits{iCond}(strcmp(trialData.session_ID,sesid));
        if sum(trialidx)>params.minNtrialscond
            
            for iNx = 1:nNeurons
                for iNy = 1:nNeurons
                    for iT = 1:params.nTimebins
                        noisecorrmat(iCond,neuroncounter + iNx,neuroncounter + iNy,iT) = corr(squeeze(tensor(iNx,iT,trialidx)),squeeze(tensor(iNy,iT,trialidx)));
                    end
                end
            end
            
        end
    end
    
    neuroncounter = neuroncounter + length(tempspikeData.session_ID);
end

%% Set all autocorrelation values of 1 to NaN
noisecorrmat(noisecorrmat>0.99)   = NaN;

%% Make the figure:
params.baselinecorr     = 1;
params.nstd             = 2;
params.idx_bsl          = params.xtime>=-.5e6 & params.xtime<=0e6;

tBin                    = NaN(2,nConditions);

for iExp = [2 3] %Only for trained task contexts:
    idx_exp             = ismember(spikeData.session_ID,sessionData.session_ID(ismember(sessionData.Experiment,params.Experiments(iExp)))); %get all neurons in this task context
    idx_exp             = idx_exp & spikeData.avg_fr>1; %get all neurons in this task context
    nCross              = sum(idx_exp)^2; %calc max combinations
    datatoplot          = NaN(nConditions,nCross,params.nTimebins); %init var
    
    for iCond = 1:nConditions %For each quantile reshape spike count correlation values:
        temp                        = noisecorrmat(iCond,idx_exp,idx_exp,:);
        temp                        = reshape(temp,sum(idx_exp)^2,params.nTimebins);
        datatoplot(iCond,:,:)       = temp;
    end
    
    idx             = any(~isnan(datatoplot),3); %calculate how many correlation values:
    fprintf('n = %d pairwise combinations\n',sum(sum(idx)))
    
    %Compute mean and sem:
    meantoplot      = squeeze(nanmean(datatoplot,2));
    errortoplot     = squeeze(nanstd(datatoplot,[],2) ./ sqrt(sum(~isnan(datatoplot),2)));
    
    if params.baselinecorr %Subtract baseline values from 
        meantoplot      = meantoplot - repmat(nanmean(meantoplot(:,params.idx_bsl),2),1,params.nTimebins);
    end
    
    %Figure of overall averaged noise correlations for the three reaction time quantiles:
    figure; hold all; set(gcf,'units','normalized','Position',[0.45 0.45 0.32 0.3],'color','w');
    handles = [];
    
    for iCond = 2:nConditions
        h = shadedErrorBar(params.xtime, meantoplot(iCond,:),errortoplot(iCond,:),{params.lines_quantiles{iCond},'LineWidth',3,'Color',params.colors_quantiles(iCond,:)},0);
        handles(end+1) = h.mainLine; %#ok<SAGROW>
        delete(h.edge(1)); delete(h.edge(2)); %h.edge(1).Color = [1 1 1];  h.edge(2).Color =  [1 1 1];
    end
        
    for iCond = 2:nConditions %Skip correct rejection, start at 2
        
        temp = NaN(1000,1);
        tempidx = find(~isnan(datatoplot(iCond,:,1)));
        if tempidx
           
            CImean = repmat(nanmean(meantoplot(iCond,params.idx_bsl)),1,params.nTimebins);
            CIband = repmat(params.nstd*nanstd(meantoplot(iCond,params.idx_bsl)),1,params.nTimebins);
           
            plot(params.xtime,CImean-CIband,'--','LineWidth',1,'Color',params.colors_quantiles(iCond,:));
            
            temp = find(meantoplot(iCond,:)<(CImean-CIband) & params.xtime>0,1);
            if temp
                tBin(iExp-1,iCond) = temp;
                plot([params.xtime(tBin(iExp-1,iCond)) params.xtime(tBin(iExp-1,iCond))],[meantoplot(iCond,tBin(iExp-1,iCond)) 0.08],'Color',params.colors_quantiles(iCond,:))
                %             plot([params.xtime(tBin(iCond)) params.xtime(tBin(iCond))],[meantoplot(iCond,tBin(iCond)) 0.08],'Color',params.colors_trialtypes{ceil(iCond/4)})
            end
        end
    end
    for iCond = 1:nConditions
        if ~isnan(tBin(iExp-1,iCond))
            fprintf('Earliest timepoint of significant decorrelation for %s, %s: %3.0f ms\n',params.ExperimentLabels{iExp},params.labels_quantiles{iCond},params.xtime(tBin(iExp-1,iCond))*1e-3)
        end
    end
    
    %Figure makeup:
    plot([0 0],[-0.2 0.2],'k','LineWidth',2)
    ylabel('Baseline corrected rsc','FontSize',20)
    legend(handles,params.labels_quantiles(2:end)); legend boxoff;
    xlim([-0.3e6 1e6])
    if params.baselinecorr
        ylim([-0.075 0.05])
    else
        ylim([-0.02 0.12])
    end
    xlabel('Time (ms)')
    set(gca,'XTick',-1.2e6:0.2e6:1e6,'XTickLabels',(-1.2e6:0.2e6:1e6)*1e-3,'YTick',[-0.1 -0.075 0 0.05 0.1])
%     grid on;
    
    %     export_fig(gcf,fullfile(params.savedir,sprintf('V1_NC_Quantiles_%s_alignStim.eps',params.ExperimentLabels{iExp})))
end

%% Correlate response latency with onset of decorrelation:

tBin(isnan(tBin)) = 1;

RT_toplot       = RT_quantile;
tDecor_toplot   = params.xtime(tBin(:,2:end));
tDecor_toplot(tDecor_toplot==params.xtime(1)) = NaN;

figure; hold all; set(gcf,'units','normalized','Position',[0.45 0.45 0.16 0.23],'color','w');
scatter(RT_toplot(1,:)*1e-3,tDecor_toplot(1,:)*1e-3,120,params.colors_quantiles(2:end,:),'o','filled')
scatter(RT_toplot(2,:)*1e-3,tDecor_toplot(2,:)*1e-3,120,params.colors_quantiles(2:end,:),'s','filled')
plot([0 1e6],[0 1e6],'k:','LineWidth',2)
xlabel('Median RT')
ylabel('Earliest Decor')
set(gca,'XTick',[0 500 1000],'YTick',[0 500 1000])
xlim([0 1000])
ylim([0 1000])

% [r,p] = corr(RT_toplot(:),tDecor_toplot(:));
idx_nan = isnan(RT_toplot) | isnan(tDecor_toplot); 
[r,p] = corr(RT_toplot(~idx_nan),tDecor_toplot(~idx_nan));
text(100,600,sprintf('r=%1.3f,p=%1.3f',r,p),'FontSize',12)
legend({'UST' 'MST'},'Location','NorthWest'); legend boxoff

%%

















%% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Show noise correlations over time for each of quantiles of reaction times but ALIGNED TO FIRST LICK
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%% Compute noise correlations:

%Create timestamp of the first lick:
trialData.firstLick         = trialData.stimChange;
trialData.firstLick(~isnan(trialData.responseLatency)) = trialData.firstLick(~isnan(trialData.responseLatency)) + trialData.responseLatency(~isnan(trialData.responseLatency));

% add pseudo reaction times for no lick trials:
z = trialData.responseLatency(trialData.responseLatency>150e3);
h = fitdist(z,'lognormal');

nNoLick                     = sum(isnan(trialData.responseLatency));
%lognrnd generates a random number from the lognormal distribution with the distribution parameters mu (mean of logarithmic values) and sigma (standard deviation of logarithmic values).
newRT                       = lognrnd(h.mu,h.sigma,[nNoLick 1]);
newRT(newRT<100e3)          = lognrnd(h.mu,h.sigma,[sum(newRT<100e3) 1]);
trialData.firstLick(isnan(trialData.responseLatency)) = trialData.stimChange(isnan(trialData.responseLatency)) + newRT;

%OLD method:
% nNoLick                     = sum(isnan(trialData.responseLatency));
% respmean                    = nanmean(trialData.responseLatency);
% respstd                     = nanstd(trialData.responseLatency);
% trialData.firstLick(isnan(trialData.responseLatency)) = trialData.stimChange(isnan(trialData.responseLatency)) + randn(nNoLick,1)*respstd+respmean;

%%
params.AlignOn              = 'firstLick';

params.t_pre                = -2.5e6;
params.t_post               = +2.5e6;

%Construct bin edges and time axis
params.edges                = [params.t_pre:params.binsize:params.t_post] - params.binsize/2;                    %#ok<NBRAK> Define edges of the bins
params.xtime                = params.t_pre:params.binsize:params.t_post-params.binsize;                        %Define time axis
params.nTimebins            = length(params.xtime); %number of time bins


noisecorrmat                = NaN(nConditions,nTotalNeurons,nTotalNeurons,params.nTimebins); %init output matrix

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
    
    %subtract mean response for each condition:
    for iCond = 1:nConditions
        trialidx                = splits{iCond}(strcmp(trialData.session_ID,sesid));
        tensor(:,:,trialidx)    = tensor(:,:,trialidx) - repmat(mean(tensor(:,:,trialidx),3),1,1,sum(trialidx));
    end
    
    %Compute noise corr:
    for iCond = 1:nConditions
        trialidx                = splits{iCond}(strcmp(trialData.session_ID,sesid));
        if sum(trialidx)>params.minNtrialscond
            
            for iNx = 1:nNeurons
                for iNy = 1:nNeurons
                    for iT = 1:params.nTimebins
                        noisecorrmat(iCond,neuroncounter + iNx,neuroncounter + iNy,iT) = corr(squeeze(tensor(iNx,iT,trialidx)),squeeze(tensor(iNy,iT,trialidx)));
                    end
                end
            end
            
        end
    end
    
    neuroncounter = neuroncounter + length(tempspikeData.session_ID);
end

%% Set all autocorrelation values of 1 to NaN
noisecorrmat(noisecorrmat>0.99)   = NaN;

%%
params.baselinecorr     = 1;
params.idx_bsl          = params.xtime>-1e6 & params.xtime<-0.3e6;
params.bsl_end          = params.xtime(find(params.idx_bsl,1,'last'));

params.nstd             = 2;

tBin                    = NaN(2,nConditions);

for iExp = [2 3] %Only for trained task contexts:
    idx_exp             = ismember(spikeData.session_ID,sessionData.session_ID(ismember(sessionData.Experiment,params.Experiments(iExp)))); %get all neurons in this task context
    idx_exp             = idx_exp & spikeData.avg_fr>1; %get all neurons in this task context
    nCross              = sum(idx_exp)^2; %calc max combinations
    datatoplot          = NaN(nConditions,nCross,params.nTimebins); %init var
    
    for iCond = 1:nConditions %For each quantile reshape spike count correlation values:
        temp                        = noisecorrmat(iCond,idx_exp,idx_exp,:);
        temp                        = reshape(temp,sum(idx_exp)^2,params.nTimebins);
        datatoplot(iCond,:,:)       = temp;
    end
    
    idx             = any(~isnan(datatoplot),3); %calculate how many correlation values:
    fprintf('n = %d pairwise combinations\n',sum(sum(idx)))
    
    %Compute mean and sem:
    meantoplot      = squeeze(nanmean(datatoplot,2));
    errortoplot     = squeeze(nanstd(datatoplot,[],2) ./ sqrt(sum(~isnan(datatoplot),2)));
    
    if params.baselinecorr %Subtract baseline values from 
        meantoplot      = meantoplot - repmat(nanmean(meantoplot(:,params.idx_bsl),2),1,params.nTimebins);
    end
    
    %Figure of overall averaged noise correlations for the three reaction time quantiles:
    figure; hold all; set(gcf,'units','normalized','Position',[0.45 0.45 0.32 0.3],'color','w');
    handles = [];
    
    for iCond = 2:nConditions
        h = shadedErrorBar(params.xtime, meantoplot(iCond,:),errortoplot(iCond,:),{params.lines_quantiles{iCond},'LineWidth',3,'Color',params.colors_quantiles(iCond,:)},0);
        handles(end+1) = h.mainLine; %#ok<SAGROW>
        delete(h.edge(1)); delete(h.edge(2)); %h.edge(1).Color = [1 1 1];  h.edge(2).Color =  [1 1 1];
    end
        
    for iCond = 2:nConditions
        
        temp        = NaN(1000,1);
        tempidx     = find(~isnan(datatoplot(iCond,:,1)));
        if tempidx
            
            CImean = repmat(nanmean(meantoplot(iCond,params.idx_bsl)),1,params.nTimebins);
            CIband = repmat(params.nstd*nanstd(meantoplot(iCond,params.idx_bsl)),1,params.nTimebins);
            
            h = plot(params.xtime,CImean-CIband,'--','LineWidth',1,'Color',params.colors_quantiles(iCond,:));
            
            temp = find(meantoplot(iCond,:)<(CImean-CIband) & params.xtime>params.bsl_end,1);
            if temp
                tBin(iExp-1,iCond) = temp;
                plot([params.xtime(tBin(iExp-1,iCond)) params.xtime(tBin(iExp-1,iCond))],[meantoplot(iCond,tBin(iExp-1,iCond)) 0.08],'Color',params.colors_quantiles(iCond,:))
            end
        end
    end
    for iCond = 1:nConditions
        if ~isnan(tBin(iExp-1,iCond))
            fprintf('Earliest timepoint of significant decorrelation for %s, %s: %3.0f ms\n',params.ExperimentLabels{iExp},params.labels_quantiles{iCond},params.xtime(tBin(iExp-1,iCond))*1e-3)
        end
    end
    
    %Figure makeup:
    plot([0 0],[-0.2 0.2],'k','LineWidth',2)
    ylabel('Baseline corrected rsc','FontSize',20)
    legend(handles,params.labels_quantiles(2:end)); legend boxoff;
    xlim([-0.6e6 1e6])
    ylim([-0.09 0.05])
    xlabel('Time (ms)')
    set(gca,'XTick',-2e6:0.2e6:2e6,'XTickLabels',(-2e6:0.2e6:2e6)*1e-3,'YTick',[-0.1 -0.05 0 0.05 0.1])
    grid on;
    
    %     export_fig(gcf,fullfile(params.savedir,sprintf('V1_NC_Quantiles_%s_alignStim.eps',params.ExperimentLabels{iExp})))
end

%%


%% Correlate response latency with onset of decorrelation:

RT_toplot       = RT_quantile;
tDecor_toplot   = params.xtime(tBin(:,2:end));

figure; hold all; set(gcf,'units','normalized','Position',[0.45 0.45 0.16 0.23],'color','w');
scatter(RT_toplot(1,:)*1e-3,tDecor_toplot(1,:)*1e-3,120,params.colors_quantiles(2:end,:),'o','filled')
scatter(RT_toplot(2,:)*1e-3,tDecor_toplot(2,:)*1e-3,120,params.colors_quantiles(2:end,:),'s','filled')
plot([0 1e6],[0 1e6],'k:','LineWidth',2)
xlabel('Median RT')
ylabel('Earliest Decor')
set(gca,'XTick',[0 500 1000],'YTick',[-500 -250 0])
xlim([0 1000])
ylim([-500 0])

[r,p] = corr(RT_toplot(:),tDecor_toplot(:));
text(100,-400,sprintf('r=%1.3f,p=%1.3f',r,p),'FontSize',12)
legend({'UST' 'MST'},'Location','NorthWest'); legend boxoff

%%
















%% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Figure supplement:
% Show noise correlations over time aligned to LICK for each of the trial types 
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%Create timestamp of the first lick:
trialData.firstLick         = trialData.stimChange;
trialData.firstLick(~isnan(trialData.responseLatency)) = trialData.firstLick(~isnan(trialData.responseLatency)) + trialData.responseLatency(~isnan(trialData.responseLatency));

% add pseudo reaction times for no lick trials:
z = trialData.responseLatency(trialData.responseLatency>150e3);
h = fitdist(z,'lognormal');

nNoLick                     = sum(isnan(trialData.responseLatency));
%lognrnd generates a random number from the lognormal distribution with the distribution parameters mu (mean of logarithmic values) and sigma (standard deviation of logarithmic values).
newRT                       = lognrnd(h.mu,h.sigma,[nNoLick 1]);
newRT(newRT<100e3)          = lognrnd(h.mu,h.sigma,[sum(newRT<100e3) 1]);
trialData.firstLick(isnan(trialData.responseLatency)) = trialData.stimChange(isnan(trialData.responseLatency)) + newRT;


%% Trial type splits:

%Parameters for trial splits:
params.labels_splits            = {'AuMiss' 'AuHit' 'VisMiss' 'VisHit' 'CR' 'FA'};
params.lines_splits             = {'-.' '-' '-.' '-' '-.' '-'};

% Get indices of the trials per trial type:
splits              = {};
splits{end+1}       = ismember(trialData.trialType,{'Y'}) & trialData.vecResponse==3;
splits{end+1}       = ismember(trialData.trialType,{'Y'}) & trialData.vecResponse==1;
splits{end+1}       = ismember(trialData.trialType,{'X'}) & trialData.vecResponse==3;
splits{end+1}       = ismember(trialData.trialType,{'X'}) & trialData.vecResponse==2;
splits{end+1}       = ismember(trialData.trialType,{'P'}) & trialData.vecResponse==3;
splits{end+1}       = ismember(trialData.trialType,{'P'}) & ismember(trialData.vecResponse,[1 2]);

%% Compute noise correlations:
%Set parameters:
params.AlignOn              = 'firstLick';

%Some measures of dimensionality:
nTotalNeurons               = length(spikeData.session_ID);
nSessions                   = length(sessionData.session_ID);

params.t_pre                = -2.5e6;
params.t_post               = +2.5e6;

%Construct bin edges and time axis
params.edges                = [params.t_pre:params.binsize:params.t_post] - params.binsize/2;                    %#ok<NBRAK> Define edges of the bins
params.xtime                = params.t_pre:params.binsize:params.t_post-params.binsize;                        %Define time axis
params.nTimebins            = length(params.xtime); %number of time bins

nConditions                 = numel(splits);
noisecorrmat                = NaN(nConditions,nTotalNeurons,nTotalNeurons,params.nTimebins); %init output matrix
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
    
    %subtract mean response for each condition:
    for iCond = 1:nConditions
        trialidx                = splits{iCond}(strcmp(trialData.session_ID,sesid));
        tensor(:,:,trialidx)    = tensor(:,:,trialidx) - repmat(mean(tensor(:,:,trialidx),3),1,1,sum(trialidx));
    end
    
    %Compute noise corr:
    for iCond = 1:nConditions
        trialidx                = splits{iCond}(strcmp(trialData.session_ID,sesid));
        if sum(trialidx)>params.minNtrialscond
            
            for iNx = 1:nNeurons
                for iNy = 1:nNeurons
                    for iT = 1:params.nTimebins
                        noisecorrmat(iCond,neuroncounter + iNx,neuroncounter + iNy,iT) = corr(squeeze(tensor(iNx,iT,trialidx)),squeeze(tensor(iNy,iT,trialidx)));
                    end
                end
            end
            
        end
    end
    
    neuroncounter = neuroncounter + length(tempspikeData.session_ID);
end

%% Set all autocorrelation values of 1 to NaN
noisecorrmat(noisecorrmat>0.99)   = NaN;

%% Make figure:
params.baselinecorr     = 0;

params.idx_bsl          = params.xtime>-1e6 & params.xtime<-0.5e6;
params.bsl_end          = params.xtime(find(params.idx_bsl,1,'last'));

params.nstd             = 2;

tBin                    = NaN(3,nConditions);

for iExp = 1:3 %For each task contingency cohort:
    idx_exp             = ismember(spikeData.session_ID,sessionData.session_ID(ismember(sessionData.Experiment,params.Experiments(iExp)))); %get all neurons in this task context
    nCross              = sum(idx_exp)^2; %calc max combinations
    datatoplot          = NaN(nConditions,nCross,params.nTimebins); %init var
    
    for iCond = 1:nConditions %For each condition reshape spike count correlation values:
        temp                        = noisecorrmat(iCond,idx_exp,idx_exp,:);
        temp                        = reshape(temp,sum(idx_exp)^2,params.nTimebins);
        datatoplot(iCond,:,:)       = temp;
    end
    
    idx             = any(~isnan(datatoplot),3); %calculate how many correlation values:
    fprintf('n = %d pairwise combinations\n',sum(sum(idx)))
    
    %Compute mean and sem:
    meantoplot      = squeeze(nanmean(datatoplot,2));
    errortoplot     = squeeze(nanstd(datatoplot,[],2) ./ sqrt(sum(~isnan(datatoplot),2)));
    
    if params.baselinecorr %Subtract baseline values 
        meantoplot      = meantoplot - repmat(nanmean(meantoplot(:,params.idx_bsl),2),1,params.nTimebins);
    end
    
    %Figure of overall averaged noise correlations for the three reaction time quantiles:
    figure; hold all; set(gcf,'units','normalized','Position',[0.45 0.45 0.32 0.3],'color','w');
    handles = [];
    
%     for iCond = 1:nConditions
    for iCond = [2 4 6]
        h = shadedErrorBar(params.xtime, meantoplot(iCond,:),errortoplot(iCond,:),{params.lines_splits{iCond},'LineWidth',3,'Color',params.colors_trialtypes{ceil(iCond/2)}},0);
        handles(end+1) = h.mainLine; %#ok<SAGROW>
        delete(h.edge(1)); delete(h.edge(2));
        %         handles(end+1) = plot(params.xtime, meantoplot(iCond,:) - nanmean(meantoplot(iCond,params.xtime<0)),params.lines_splits{iCond},'LineWidth',4,'Color',params.colors_splits(iCond,:));
    end
    
    for iCond = [2 4 6]
        
        temp        = NaN(1000,1);
        tempidx     = find(~isnan(datatoplot(iCond,:,1)));
        if tempidx
            
            CImean = repmat(nanmean(meantoplot(iCond,params.idx_bsl)),1,params.nTimebins);
            CIband = repmat(params.nstd*nanstd(meantoplot(iCond,params.idx_bsl)),1,params.nTimebins);
            
            h = plot(params.xtime,CImean-CIband,'--','LineWidth',1,'Color',params.colors_trialtypes{ceil(iCond/2)});
            
            temp = find(meantoplot(iCond,:)<(CImean-CIband) & params.xtime>params.bsl_end,1);
            if temp
                tBin(iExp,iCond) = temp;
                plot([params.xtime(tBin(iExp,iCond)) params.xtime(tBin(iExp,iCond))],[meantoplot(iCond,tBin(iExp,iCond)) 0.13],'Color',params.colors_trialtypes{ceil(iCond/2)})
            end
        end
    end
    for iCond = 1:nConditions
        if ~isnan(tBin(iExp,iCond))
            fprintf('Earliest timepoint of significant decorrelation for %s, %s: %3.0f ms\n',params.ExperimentLabels{iExp},params.labels_splits{iCond},params.xtime(tBin(iExp,iCond))*1e-3)
        end
    end
    
    %Figure makeup:
    plot([0 0],[-0.2 0.2],'k','LineWidth',2)
    ylabel('Baseline corrected rsc','FontSize',20)
    legend(handles,params.labels_splits([2 4 6])); legend boxoff;
%     xlim([-0.2e6 1.2e6])
    ylim([0 0.13])
    
    xlim([-0.5e6 1e6])
    xlabel('Time (ms)')
    set(gca,'XTick',-2e6:0.2e6:2e6,'XTickLabels',(-2e6:0.2e6:2e6)*1e-3,'YTick',[-0.1 -0.05 0 0.05 0.1])
    set(gca,'XTick',-2.1e6:0.3e6:2.1e6,'XTickLabels',(-2.1e6:0.3e6:2.1e6)*1e-3,'YTick',[-0.1 -0.05 0 0.05 0.1])
    grid on;
    
    if params.baselinecorr %Subtract baseline values from
        ylim([-0.07 0.05])
    end
    export_fig(gcf,fullfile(params.savedir,sprintf('V1_NC_trialtypes_%s_alignLick.eps',params.ExperimentLabels{iExp})))
end




