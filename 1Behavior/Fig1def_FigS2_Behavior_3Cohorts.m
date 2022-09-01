%% Script that analyzes primary behavioral measures of performance across the three task versions:
% This script analyzes the behavior of animals in the audiovisual change detection task across the 3 cohorts (NE, UST, MST)
% E.g. it plots the session-averaged psychometric performance of each mouse split by cohort
% analyzes dprime, reaction time and thesholds
% Oude Lohuis et al. 2022 Nat Comms

%% Parameters
params.Experiments          = {'ChangeDetectionConflictDecor' 'VisOnlyPsychophysics' {'BehaviorConflict' 'BehaviorPsychophysics'}};
params.ExperimentLabels     = {'NE' 'UST' 'MST'};
params.nExperiments         = length(params.Experiments);

params                      = MOL_getColors_CHDET(params);

params.savedir              = 'E:\Documents\PhD\Figures\Project CHDET\Results - 2nd Bump\1 Behavior\Supp Fig 1 - Behavior Hz Oct\';

%% Load the dataset:
load('E:\Matlab\MOL_Analysis\AOudeLohuisetal_2022_NatComms\1Behavior\Data1_2.mat','sessionData','trialData')
fprintf('Dataset: %d sessions, %d trials\n',length(sessionData.session_ID),length(trialData.session_ID));

%% Fit each session:
nSessions               = length(sessionData.mousename);

FullParam               = NaN(nSessions,8); %Init output variable

for iSes = 1:nSessions %Fit each session
    fprintf('\nFitting session %d/%d \n\n',iSes,nSessions);
    %Get the data for this session only:
    [tempsessionData,temptrialData]                 = MOL_getTempPerSes(sessionData.session_ID(iSes),sessionData,trialData);
    %Get response rates per condition:
    [visconditions,auconditions,FullRespMat,~]      = MOL_Psy_GetTrialMat(tempsessionData,temptrialData);
    %Correct dimensions for some sessions:
    if numel(visconditions)==5 && numel(auconditions)==4
        FullRespMat = [FullRespMat; FullRespMat(end,:,:)];
        auconditions = [2/256 8/256 32/256 64/256 128/256];
    end
    if numel(visconditions)==6 && numel(auconditions)==4
        FullRespMat = [FullRespMat; FullRespMat(end-1:end,:,:)];
        auconditions = [1/256 2/256 8/256 32/256 64/256 128/256];
    end
    if numel(visconditions)==4 && numel(auconditions)==3
        FullRespMat = [FullRespMat(1:2,:,:); FullRespMat(2:end,:,:)];
        auconditions = [2/256 8/256 32/256 128/256];
    end
    if numel(visconditions)==5 && numel(auconditions)==8
        idx = [1 3 5 7 8];
        FullRespMat = FullRespMat([1 idx+1],:,:);
        auconditions = auconditions(idx);
    end
    
    %Rewrite fullresponse table to contingency table for model fitting (without conflict trials)
    ctable              = NaN(3,3,numel(visconditions));
    ctable(1,:,:)       = permute(FullRespMat(2:end,1,:),[3 1 2]);              %Auditory
    ctable(2,:,:)       = permute(FullRespMat(1,2:end,:),[1 3 2]);              %Visual
    ctable(3,:,:)       = repmat(permute(FullRespMat(1,1,:),[1 3 2]),1,1,numel(visconditions));    %Probe trials
    
    %Align visual and auditory intensities
    %visual and auditory conditions should be normalized such that
    %value of 1 corresponds to expected asymptotic dmax
    normconditions = [auconditions / max(auconditions); visconditions / max(visconditions)];
    
    %Fit psychometric 2ADC model:
%     fprintf('Fitting session %d/%d, of animal %d/%d\n',ses,length(sesselec),mou,length(mouseids));
    [theta_est, theta_err, LLF, ctable_mod, ctable_fit, sivals, psyc_mate, ce] = mADC_model_fit_psyc_editMOL(ctable,normconditions,[]);
    
    %Store parameters:
    %theta_est = 8 parameters: % 3 d' -related parameters and one 'c' parameter for each location
    % order is daud, dvis, naud, nvis, s50aud, s50vis, caud, vis
    % dmax = theta(1:M); n = theta(M+1:2*M); s50 = theta(2*M+1:3*M); c = theta(3*M+1:end);
    theta_est(5) = theta_est(5) * max(auconditions); %undo normalization of conditions
    theta_est(6) = theta_est(6) * max(visconditions); %undo normalization of conditions
    FullParam(iSes,:)            = theta_est; %store parameters

end

%% Report mean number of trials:
fprintf('\nAverage number of trials: %4.2f\n',nanmean(sessionData.totalTrials))
fprintf('Minimum number of trials: %4.0f\n',min(sessionData.totalTrials))
fprintf('Maximum number of trials: %4.0f\n\n',max(sessionData.totalTrials))

%% Plot average rates for each cohort:
for iExp = 1:params.nExperiments
    
    figure; set(gcf,'color','w','units','normalized','Position', [0.1 0.1 .7 .7]); hold all;
    
    expanimals                 = unique(sessionData.mousename(ismember(sessionData.Experiment,params.Experiments{iExp})));
    
    for iAnimal = 1:length(expanimals)
        
        [tempsessionData,~]        = MOL_getTempPerSes(sessionData.session_ID(strcmp(sessionData.mousename,expanimals{iAnimal})),sessionData,trialData);
        
        idx_ses  = strcmp(sessionData.mousename,expanimals{iAnimal});
        meantheta_est             = median(FullParam(idx_ses,:),1);
        
        % Settings:
        switch tempsessionData.auChangeUnit{1}
            case 'Hz'
                params.auprobepos = 0.5;
                params.auticks = [100 4000];
                params.auticklabels = ['Probe' num2cell(params.auticks)];
                params.auxaxislabel  = 'Delta frequency (Hz)';
                params.auystatslabel = 'Auditory threshold (Hz)';
            case 'Oct'
                params.auprobepos = 0.001;
                params.auticks = [1/256 1/64 1/8 1/2];
                params.auticklabels = {'Probe' '1/256' '1/64' '1/8' '1/2'};
                params.auxaxislabel  = 'Delta frequency (Oct)';
                params.auystatslabel = 'Auditory threshold (partial octave)';
        end
        
        params.visprobepos     = 0.5;
        params.visticks        = [5 30 90];
        params.vistickslabels  = ['Probe' num2cell(params.visticks)];
        params.visxaxislabel   = 'Delta orientation (Degrees)';
        params.visystatslabel  = 'Visual threshold (Degrees)';
        
        params.yticks          = [0 0.25 0.5 0.75 1];
        
        % Generate contingency table from fitted parameters:
        [~,xvals_fit_vis,ctable_fit_mat] = MOL_Gen2ADC_PsyCurve(meantheta_est,params);
        %Make common x-axis for the auditory to merge from both cohorts:
        xvals_fit_au = 10.^(linspace(log10(0.001),log10(1/2),1000)); %logarithmic spacing
        
        %Audio:
        subplot(1,2,1); hold all;
        plot(xvals_fit_au,squeeze(ctable_fit_mat(1,1,:)),'-','Color',[1 0 0],'LineWidth',3);
        plot(xvals_fit_au,squeeze(ctable_fit_mat(1,2,:)),':','Color',[0.2 0.2 1],'LineWidth',3);
        
        %Visual:
        subplot(1,2,2); hold all;
        plot(xvals_fit_vis,squeeze(ctable_fit_mat(2,2,:)),'-','Color',[0 0 1],'LineWidth',3);
        plot(xvals_fit_vis,squeeze(ctable_fit_mat(2,1,:)),':','Color',[1 0.2 0.2],'LineWidth',3);
    end
    
    params.auprobepos = 0.001;
    params.auticks = [1/256 1/64 1/8 1/2];
    params.auticklabels = {'Probe' '1/256' '1/64' '1/8' '1/2'};
    params.auxaxislabel  = 'Delta frequency (Oct)';
    params.auystatslabel = 'Auditory threshold (partial octave)';
    
    MOL_Psy2Sided_FigMakeup(params)
    
end

%% Show Dprime for all experiments
%Get all the dprimes in one multidimensional variable:
datatoplot      = NaN(3,2,20,20); %init output var (dim = experiments, modality, mouse, session)
for iExp = 1:3
    expanimals                 = unique(sessionData.mousename(ismember(sessionData.Experiment,params.Experiments{iExp})));
    
    for iAnimal = 1:length(expanimals)
        sesidx                  = ismember(sessionData.mousename,expanimals(iAnimal));
        datatoplot(iExp,1,iAnimal,1:sum(sesidx))        = FullParam(sesidx,1);
        datatoplot(iExp,2,iAnimal,1:sum(sesidx))        = FullParam(sesidx,2);
    end
end

datatoplot_re       = reshape(datatoplot,3,2,20*20); %reshape
mediandatatoplot    = nanmedian(datatoplot,4); %average over sessions

figure; set(gcf,'color','w','units','normalized','Position', [0.1 0.4 .2 .4]); hold all;
params.offset       = 1; %Offset of bars relative to each other

for iExp= 1:3
    for iMod = 1:2
        %With std across all sessions:
        %         errorbar(1+params.offset*(iExp-1)+params.offset/3*(iMod-1),nanmean(meandatatoplot(iExp,iMod,:),3),nanstd(meandatatoplot(iExp,iMod,:),[],3),'k.','MarkerSize',1,'LineWidth',5)
        %         errorbar(1+params.offset*(iExp-1)+params.offset/3*(iMod-1),nanmean(meandatatoplot(iExp,iMod,:),3),nanstd(meandatatoplot(iExp,iMod,:),[],3),'k.','MarkerSize',1,'LineWidth',5)
        %with interquartile range across animals:
        prc25       = prctile(mediandatatoplot(iExp,iMod,:),25,3) - nanmedian(mediandatatoplot(iExp,iMod,:),3);
        prc75       = prctile(mediandatatoplot(iExp,iMod,:),75,3) - nanmedian(mediandatatoplot(iExp,iMod,:),3);
        errorbar(1+params.offset*(iExp-1)+params.offset/3*(iMod-1),nanmedian(mediandatatoplot(iExp,iMod,:),3),prc25,prc75,'k.','MarkerSize',1,'LineWidth',5)
        errorbar(1+params.offset*(iExp-1)+params.offset/3*(iMod-1),nanmedian(mediandatatoplot(iExp,iMod,:),3),0,'k.','MarkerSize',40,'LineWidth',5)
    end
end

for iExp= 1:3
    for iMod = 1:2
%         scatter(1+params.offset*(iExp-1)+params.offset/3*(iMod-1)+0.15 + randn(sum(~isnan(meandatatoplot(iExp,iMod,:))),1)*0.02,meandatatoplot(iExp,iMod,~isnan(meandatatoplot(iExp,iMod,:))),50,'k','filled')
%         scatter(1+params.offset*(iExp-1)+params.offset/3*(iMod-1)+0.15 + randn(sum(~isnan(mediandatatoplot(iExp,iMod,:))),1)*0.02,mediandatatoplot(iExp,iMod,~isnan(mediandatatoplot(iExp,iMod,:))),50,params.colors_modalities{iMod},'filled')
        scatter(1+params.offset*(iExp-1)+params.offset/3*(iMod-1)+0.15 + randn(sum(~isnan(mediandatatoplot(iExp,iMod,:))),1)*0.02,...
            squeeze(mediandatatoplot(iExp,iMod,~isnan(mediandatatoplot(iExp,iMod,:)))),50,params.colors_modalities{iMod},'filled')
    end
end

ylabel('Dprime')
set(gca,'XTick', 1:3,'XTickLabels',params.ExperimentLabels,'XTickLabelRotation',45)
ylim([-0.2 3])
xlim([0.5 3.5])

% %% OLD Statistics: (Comparing visual d-prime for unisensory vs multisensory)
% Xdata = datatoplot_re(2,2,:); Xdata = Xdata(:);
% Ydata = datatoplot_re(3,2,:); Ydata = Ydata(:);
% prank = ranksum(Xdata,Ydata) * 6;
% fprintf('UST vs MST visual dprime: Rank-sum test (n=%d UST sessions, n=%d, MST sessions): p=%1.2e\n',sum(~isnan(Xdata)),sum(~isnan(Ydata)),prank)
% sigstar([2.33 3.33],prank)
% 
%Statistics: (Comparing visual d-prime with auditory dprime within multisensory trained mice)
% Xdata = datatoplot_re(3,1,:); Xdata = Xdata(:);
% Ydata = datatoplot_re(3,2,:); Ydata = Ydata(:);
% prank = signrank(Xdata,Ydata) * 3;
% fprintf('Visual vs audio dprime: wilcoxon signed rank test of different dprime (n=%d sessions): p=%1.2e\n',sum(~isnan(Ydata)),prank)
% sigstar([3 3.33],prank)

%% Multi-level statistics: 
X_coh       = NaN(nSessions,1);
X_coh(ismember(sessionData.Experiment,params.Experiments{1}))           = 1;
X_coh(ismember(sessionData.Experiment,params.Experiments{2}))           = 2;
X_coh(ismember(sessionData.Experiment,params.Experiments{3}))           = 3;

G_mou       = cell(nSessions,1);
uMice       = unique(sessionData.mousename);
for iMouse = 1:length(uMice)
%     G_mou(strcmp(sessionData.mousename,uMice{iMouse})) = iMouse;
    G_mou(strcmp(sessionData.mousename,uMice{iMouse})) = uMice(iMouse);
end

%% Comparing visual d-prime for unisensory vs multisensory
Y               = FullParam(:,2); %test visual dprime
idx             = ismember(X_coh,[2 3]); %only compare UST and MST
tbl             = table(Y(idx),X_coh(idx),G_mou(idx),'VariableNames',{'Dprime','Cohort','Mouse'}); %Create table for mixed model
lme             = fitlme(tbl,'Dprime~Cohort+(1|Mouse)'); %construct linear mixed effects model with fixed effect of cohort and random intercept for different mice
stats           = dataset2table(anova(lme,'DFMethod','Satterthwaite')); %Perform ANOVA on model and output as matrix
fprintf('Visual dprime UST vs MST : (Linear Mixed Model)\n')
fprintf('(F(%d,%2.0f) = %1.2f, p=%1.2f, ANOVA)\n',stats{2,3},stats{2,4},stats{2,2},stats{2,5})
sigstar([2.33 3.33],stats{2,5})

tempfile = fullfile(params.savedir,'SourceData_Fig1f_VisDprime_CohortUSTMST.xlsx');
writetable(tbl,tempfile)

%% Comparing visual d-prime with auditory dprime within multisensory trained mice
Y               = [FullParam(:,1); FullParam(:,2)]; %get data of visual and audio dprime
X_mod           = [ones(nSessions,1); ones(nSessions,1)*2]; %construct predictor variable of modality
G_mou_2         = [G_mou; G_mou]; %double to get index of only MST mice

idx             = [ismember(X_coh,3); ismember(X_coh,3)]; %index of only MST

tbl             = table(Y(idx),X_mod(idx),G_mou_2(idx),'VariableNames',{'Dprime','Modality','Mouse'}); %Create table for mixed model
lme             = fitlme(tbl,'Dprime~Modality+(1|Mouse)'); %construct linear mixed effects model with fixed effect of modality and random intercept for different mice
stats           = dataset2table(anova(lme,'DFMethod','Satterthwaite')); %Perform ANOVA on model and output as matrix
fprintf('Visual vs audio dprime: (Linear Mixed Model)\n')
fprintf('(F(%d,%2.0f) = %1.2f, p=%1.2e, ANOVA)\n',stats{2,3},stats{2,4},stats{2,2},stats{2,5})
sigstar([3 3.33],stats{2,5})

tempfile = fullfile(params.savedir,'SourceData_Fig1f_Dprime_Modality_MST.xlsx');
writetable(tbl,tempfile)

%% Compare thresholds visonly and multisensory:
%Get the visual threshold for trained animals in one multidimensional variable:
datatoplot = NaN(2,20,20); %init output var (dim = experiments, mouse, session)
for iExp = 2:3 %only visual and multisensory trained animals
    expanimals                 = unique(sessionData.mousename(ismember(sessionData.Experiment,params.Experiments{iExp})));
    
    for iAnimal = 1:length(expanimals)
        sesidx                  = ismember(sessionData.mousename,expanimals(iAnimal));
        datatoplot(iExp-1,iAnimal,1:sum(sesidx)) = FullParam(sesidx,6);
    end
end

%Discard sessions with absurd fitted values, e.g. threshold below 1 degree of orientation change
fprintf('%d sessions excluded due to low threshold (<1deg)\n',sum(sum(sum(datatoplot<1))));
fprintf('%d sessions excluded due to high threshold (>45deg)\n',sum(sum(sum(datatoplot>45))));
datatoplot(datatoplot<1) = NaN; datatoplot(datatoplot>45) = NaN;

% Make figure:
figure; set(gcf,'color','w','units','normalized','Position', [0.1 0.4 .2 .4]); hold all;

datatoplot_re       = reshape(datatoplot,2,20*20); %reshape all session below each other
mediandatatoplot    = nanmedian(datatoplot,3); %average over sessions

%with interquartile range and median across animals:
prc25               = prctile(datatoplot_re,25,2);
prc75               = prctile(datatoplot_re,75,2);
errorbar([1 1+params.offset],nanmedian(datatoplot_re,2),prc25-nanmedian(datatoplot_re,2),prc75-nanmedian(datatoplot_re,2),'k.','MarkerSize',1,'LineWidth',5)
errorbar([1 1+params.offset],nanmedian(datatoplot_re,2),[0 0],'k.','MarkerSize',40,'LineWidth',5)
%Plot individual animals:
scatter(1-0.2 + randn(sum(~isnan(mediandatatoplot(1,:))),1)*0.05,mediandatatoplot(1,~isnan(mediandatatoplot(1,:))),50,params.colors_experiments{2},'filled')
scatter(1-0.2+params.offset + randn(sum(~isnan(mediandatatoplot(2,:))),1)*0.05,mediandatatoplot(2,~isnan(mediandatatoplot(2,:))),50,params.colors_experiments{3},'filled')

ylabel('Visual Threshold (deg)')
set(gca,'XTick', [0.8 0.8+params.offset],'XTickLabels',params.ExperimentLabels(2:3),'XTickLabelRotation',45)


%% Statistics: (Comparing visual threshold between UST and MST trained mice)
% Xdata = datatoplot(1,:,:); Xdata = Xdata(:);
% Ydata = datatoplot(2,:,:); Ydata = Ydata(:);
% 
% prank = ranksum(Xdata,Ydata);
% 
% fprintf('Threshold: Rank-sum test of different thresholds (n=%d unisensory sessions, n=%d, multisensory sessions): p=%1.2f\n',sum(~isnan(Xdata)),sum(~isnan(Ydata)),prank)
% sigstar([1 2],prank)

%% Multi-level statistics: 
Y               = FullParam(:,6); %test visual threshold
idx             = ismember(X_coh,[2 3]); %only compare UST and MST
tbl             = table(Y(idx),X_coh(idx),G_mou(idx),'VariableNames',{'Threshold','Cohort','Mouse'}); %Create table for mixed model
lme             = fitlme(tbl,'Threshold~Cohort+(1|Mouse)'); %construct linear mixed effects model with fixed effect of cohort and random intercept for different mice
stats           = dataset2table(anova(lme,'DFMethod','Satterthwaite')); %Perform ANOVA on model and output as matrix

fprintf('Visual threshold UST vs MST : (Linear Mixed Model ANOVA)\n')
fprintf('(F(%d,%2.0f)=%1.2f, p=%1.2f, ANOVA)\n',stats{2,3},stats{2,4},stats{2,2},stats{2,5})
sigstar([1 2],stats{2,5})

tempfile = fullfile(params.savedir,'SourceData_Fig1g_VisThr_CohortUSTMST.xlsx');
writetable(tbl,tempfile)

%% Compare sensitivity visonly and multisensory:
%Get the visual threshold for trained animals in one multidimensional variable:
datatoplot = NaN(2,20,20); %init output var (dim = experiments, mouse, session)
for iExp = 2:3 %only visual and multisensory trained animals
    expanimals                 = unique(sessionData.mousename(ismember(sessionData.Experiment,params.Experiments{iExp})));
    
    for iAnimal = 1:length(expanimals)
        sesidx                  = ismember(sessionData.mousename,expanimals(iAnimal));
        datatoplot(iExp-1,iAnimal,1:sum(sesidx)) = FullParam(sesidx,4);
    end
end

% Make figure:
figure; set(gcf,'color','w','units','normalized','Position', [0.1 0.4 .2 .4]); hold all;

datatoplot_re           = reshape(datatoplot,2,20*20); %reshape all session below each other
mediandatatoplot        = nanmedian(datatoplot,3); %average over sessions

%with interquartile range and median across animals:
prc25               = prctile(datatoplot_re,25,2);
prc75               = prctile(datatoplot_re,75,2);
errorbar([1 1+params.offset],nanmedian(datatoplot_re,2),prc25-nanmedian(datatoplot_re,2),prc75-nanmedian(datatoplot_re,2),'k.','MarkerSize',1,'LineWidth',5)
errorbar([1 1+params.offset],nanmedian(datatoplot_re,2),[0 0],'k.','MarkerSize',40,'LineWidth',5)
%Plot individual animals:
scatter(1-0.2 + randn(sum(~isnan(mediandatatoplot(1,:))),1)*0.05,mediandatatoplot(1,~isnan(mediandatatoplot(1,:))),50,params.colors_experiments{2},'filled')
scatter(1-0.2+params.offset + randn(sum(~isnan(mediandatatoplot(2,:))),1)*0.05,mediandatatoplot(2,~isnan(mediandatatoplot(2,:))),50,params.colors_experiments{3},'filled')

ylabel('Sensitivity (n)')
set(gca,'XTick', [0.8 0.8+params.offset],'XTickLabels',params.ExperimentLabels(2:3),'XTickLabelRotation',45)
set(gca,'YScale','log')



%% OLD Statistics: (Comparing visual threshold between UST and MST trained mice)

% Xdata = datatoplot(1,:,:); Xdata = Xdata(:);
% Ydata = datatoplot(2,:,:); Ydata = Ydata(:);
% prank = ranksum(Xdata,Ydata);
% % fprintf('Sensitivity: Rank-sum test of different thresholds (n=%d unisensory sessions, n=%d, multisensory sessions): p=%1.2f\n',sum(~isnan(Xdata)),sum(~isnan(Ydata)),prank)
% sigstar([1 2],prank)

%% NEW Statistics: Multi-level statistics: 
Y               = FullParam(:,4); %test visual threshold
idx             = ismember(X_coh,[2 3]); %only compare UST and MST
tbl             = table(Y(idx),X_coh(idx),G_mou(idx),'VariableNames',{'Threshold','Cohort','Mouse'}); %Create table for mixed model
lme             = fitlme(tbl,'Threshold~Cohort+(1|Mouse)'); %construct linear mixed effects model with fixed effect of cohort and random intercept for different mice
stats           = dataset2table(anova(lme,'DFMethod','Satterthwaite')); %Perform ANOVA on model and output as matrix

fprintf('Visual sensitivity UST vs MST : (Linear Mixed Model ANOVA)\n')
fprintf('(F(%d,%2.0f) = %1.2f, p = %1.2f, ANOVA)\n',stats{2,3},stats{2,4},stats{2,2},stats{2,5})
sigstar([1 2],stats{2,5})

tempfile = fullfile(params.savedir,'SourceData_Fig1g_VisSens_CohortUSTMST.xlsx');
writetable(tbl,tempfile)

%% Make figure of median response latency:
splits              = {};
params.trialcolors  = {};

%Take only responses after 200ms (grace period) and in the same experimental setup (to control for any differences between boxes)
idx_lat             = ~(strcmp(trialData.trialType,'X') & trialData.responseLatency<200e3) | isnan(trialData.responseLatency);
idx_box             = ismember(trialData.session_ID,unique(sessionData.session_ID(sessionData.Box==9)));

%UST visual only:
iExp                = 2; %select only trials from sessions in UST mice:
idx_exp             = ismember(trialData.session_ID,sessionData.session_ID(ismember(sessionData.Experiment,params.Experiments{iExp})));

for iSplit = 2:5%visual trials of different intensities:
    splits{end+1}               = strcmp(trialData.trialType,'X') & trialData.visualOriChangeNorm==iSplit & trialData.vecResponse==2 & idx_exp & idx_box & idx_lat;
    params.trialcolors{end+1}   = [0 0 iSplit/5];
end

%MST visual and auditory:
iExp                = 3; %select only trials from sessions in MST mice:
idx_exp             = ismember(trialData.session_ID,sessionData.session_ID(ismember(sessionData.Experiment,params.Experiments{iExp})));

for iSplit = 2:5 %visual trials of different intensities:
    splits{end+1}               = strcmp(trialData.trialType,'X') & trialData.visualOriChangeNorm==iSplit & trialData.vecResponse==2 & idx_exp & idx_box & idx_lat;
    params.trialcolors{end+1}   = [0 0 iSplit/5];
end

for iSplit = 2:5%audio trials of different intensities:
    splits{end+1}               = strcmp(trialData.trialType,'Y') & trialData.audioFreqChangeNorm==iSplit & trialData.vecResponse==1 & idx_exp & idx_box & idx_lat; %#ok<*AGROW>
    params.trialcolors{end+1}   = [iSplit/5 0 0];
end

params.triallabels = {'Vsub' 'Vthr' 'Vsup' 'Vmax'  'Vsub' 'Vthr' 'Vsup' 'Vmax' 'Asub' 'Athr' 'Asup' 'Amax' };
%Init output vars:
nSplits                 = length(splits);
datatoplot              = NaN(nSplits,10000);
mediantoplot            = NaN(nSplits,1);
lowererror              = NaN(nSplits,1);
uppererror              = NaN(nSplits,1);
% errortoplot             = NaN(nSplits,1);

%For each trial type store all the response latencies and compute mean and sem:
for iSplit = 1:nSplits
    datatoplot(iSplit,1:sum(splits{iSplit})) = trialData.responseLatency(splits{iSplit});
    mediantoplot(iSplit)    = nanmedian(datatoplot(iSplit,:));
    lowererror(iSplit)      = prctile(datatoplot(iSplit,:),25) - nanmedian(datatoplot(iSplit,:));
    uppererror(iSplit)      = prctile(datatoplot(iSplit,:),75) - nanmedian(datatoplot(iSplit,:));
%      lowererror(iSplit)   = mad(datatoplot(iSplit,:));
%     errortoplot(iSplit)   = nanstd(datatoplot(iSplit,:));
%     errortoplot(iSplit)   = errortoplot(iSplit) / sqrt(sum(splits{iSplit}));
end

mediantoplot            = mediantoplot*1e-3; %convert to ms
lowererror              = lowererror*1e-3;
uppererror              = uppererror*1e-3;
% errortoplot             = errortoplot*1e-3;

figure; hold all; set(gcf,'units','normalized','Position',[0.1 0.25 0.18 0.28],'color','w');
h = boxplot(datatoplot'*1e-3,'plotstyle','compact','whisker',0.5,'outliersize',0.0001,'widths',10);
for iSplit = 1:length(splits)
    set(h(1,iSplit),'Color',params.trialcolors{iSplit})
    set(h(2,iSplit),'Color',params.trialcolors{iSplit})
    set(h(3,iSplit),'MarkerEdgeColor',params.trialcolors{iSplit})
    set(h(4,iSplit),'MarkerEdgeColor',params.trialcolors{iSplit})
%     h(iSplit)   = bar(iSplit,mediantoplot(iSplit),0.8);
%     set(h(iSplit),'facecolor',params.trialcolors{iSplit});
%     errorbar(iSplit,mediantoplot(iSplit),lowererror(iSplit),uppererror(iSplit),'k.')
%     boxplot(iSplit,datatoplot(iSplit,:),'plotstyle','compact')
%     errorbar(iSplit,mediantoplot(iSplit),errortoplot(iSplit),'k.')
end
% h = violi(datatoplot'*1e-3,'plotstyle','compact');

set(gca,'XTick',1:nSplits,'XTickLabel',params.triallabels,'XTickLabelRotation',45)

xlabel('')
ylabel('Time (ms)')
xlim([0 nSplits+1])
ylim([100 1000])

%% OLD Statistical testing:
%Reorder data:
% nGroups             = 12; %12 conditions (3 modalities -vision UST, vision and audition MST- 4 per modality)
% nInit               = 10000;
% groups              = repmat(1:nGroups,nInit,1); groups = reshape(groups,nGroups*nInit,1);
% datatotest          = reshape(datatoplot',nGroups*nInit,1); %reshape to one column vector
% groups              = groups(~isnan(datatotest)); %filter out nans
% datatotest          = datatotest(~isnan(datatotest)); %filter out nans

%% NEW Multilevel statistical testing:

Y_rt                = trialData.responseLatency;
X_sal_vis           = trialData.visualOriChangeNorm-1;

nTrials             = length(trialData.session_ID);
G_mou               = cell(nTrials,1);

for iMouse = 1:length(uMice)
%     G_mou(ismember(trialData.session_ID,sessionData.session_ID(strcmp(sessionData.mousename,uMice{iMouse})))) = iMouse;
    G_mou(ismember(trialData.session_ID,sessionData.session_ID(strcmp(sessionData.mousename,uMice{iMouse})))) = uMice(iMouse);
end

%Visual saliency on RT (UST)
iExp                = 2;
idx_exp             = ismember(trialData.session_ID,sessionData.session_ID(ismember(sessionData.Experiment,params.Experiments{iExp})));

idx                 = strcmp(trialData.trialType,'X') & ismember(trialData.visualOriChangeNorm,2:5) & trialData.vecResponse==2 & idx_exp & idx_box & idx_lat;

tbl                 = table(log(Y_rt(idx)),X_sal_vis(idx),G_mou(idx),'VariableNames',{'RT','Saliency','Mouse'}); %Create table for mixed model
lme                 = fitlme(tbl,'RT~Saliency+(1|Mouse)'); %construct linear mixed effects model with fixed effect of cohort and random intercept for different mice
stats               = dataset2table(anova(lme,'DFMethod','Satterthwaite')); %Perform ANOVA on model and output as matrix

fprintf('Fixed effect of visual saliency on RT (UST) (%d trials): (Linear Mixed Model ANOVA)\n',sum(idx))
fprintf('(F(%d,%2.0f) = %1.2f, p=%1.2f, ANOVA)\n',stats{2,3},stats{2,4},stats{2,2},stats{2,5})

tempfile = fullfile(params.savedir,'SourceData_Fig1h_VisSal_RT_UST.xlsx');
writetable(tbl,tempfile)

%Visual saliency on RT (MST)
iExp                = 3;
idx_exp             = ismember(trialData.session_ID,sessionData.session_ID(ismember(sessionData.Experiment,params.Experiments{iExp})));

idx                 = strcmp(trialData.trialType,'X') & ismember(trialData.visualOriChangeNorm,2:5) & trialData.vecResponse==2 & idx_exp & idx_box & idx_lat;

tbl                 = table(Y_rt(idx),X_sal_vis(idx),G_mou(idx),'VariableNames',{'RT','Saliency','Mouse'}); %Create table for mixed model
% tbl                 = table(log(Y_rt(idx)),X_sal(idx),G_mou(idx),'VariableNames',{'RT','Saliency','Mouse'}); %Create table for mixed model
lme                 = fitlme(tbl,'RT~Saliency+(1|Mouse)'); %construct linear mixed effects model with fixed effect of cohort and random intercept for different mice
stats               = dataset2table(anova(lme,'DFMethod','Satterthwaite')); %Perform ANOVA on model and output as matrix

fprintf('Fixed effect of visual saliency on RT (MST) (%d trials): (Linear Mixed Model ANOVA)\n',sum(idx))
fprintf('(F(%d,%2.0f) = %1.2f, p=%1.2e, ANOVA)\n',stats{2,3},stats{2,4},stats{2,2},stats{2,5})

tempfile = fullfile(params.savedir,'SourceData_Fig1h_VisSal_RT_MST.xlsx');
writetable(tbl,tempfile)

%Auditory saliency on RT (MST)
iExp                = 3;
idx_exp             = ismember(trialData.session_ID,sessionData.session_ID(ismember(sessionData.Experiment,params.Experiments{iExp})));
X_sal_aud           = trialData.audioFreqChangeNorm-1;

idx                 = strcmp(trialData.trialType,'Y') & ismember(trialData.audioFreqChangeNorm,2:5) & trialData.vecResponse==1 & idx_exp & idx_box & idx_lat;

tbl                 = table(Y_rt(idx),X_sal_aud(idx),G_mou(idx),'VariableNames',{'RT','Saliency','Mouse'}); %Create table for mixed model
% tbl                 = table(log(Y_rt(idx)),X_sal(idx),G_mou(idx),'VariableNames',{'RT','Saliency','Mouse'}); %Create table for mixed model
lme                 = fitlme(tbl,'RT~Saliency+(1|Mouse)'); %construct linear mixed effects model with fixed effect of cohort and random intercept for different mice
stats               = dataset2table(anova(lme,'DFMethod','Satterthwaite')); %Perform ANOVA on model and output as matrix

fprintf('Fixed effect of auditory saliency on RT (MST) (%d trials): (Linear Mixed Model ANOVA)\n',sum(idx))
fprintf('(F(%d,%2.0f) = %1.2f, p=%1.2e, ANOVA)\n',stats{2,3},stats{2,4},stats{2,2},stats{2,5})

tempfile = fullfile(params.savedir,'SourceData_Fig1h_AuSal_RT_MST.xlsx');
writetable(tbl,tempfile)

%Modality on RT (MST):
iExp                = 3;
idx_exp             = ismember(trialData.session_ID,sessionData.session_ID(ismember(sessionData.Experiment,params.Experiments{iExp})));

X_sal               = X_sal_aud + X_sal_vis; %combine vectors of saliency
idx_vis             = strcmp(trialData.trialType,'X') & ismember(trialData.visualOriChangeNorm,2:5) & trialData.vecResponse==2 & idx_exp & idx_box & idx_lat;
idx_aud             = strcmp(trialData.trialType,'Y') & ismember(trialData.audioFreqChangeNorm,2:5) & trialData.vecResponse==1 & idx_exp & idx_box & idx_lat;

idx                 = idx_vis | idx_aud;
X_mod               = zeros(size(idx)); X_mod(idx_vis)=1; X_mod(idx_aud)=2;

tbl                 = table(Y_rt(idx),X_mod(idx),X_sal(idx),G_mou(idx),'VariableNames',{'RT','Modality','Saliency','Mouse'}); %Create table for mixed model
% tbl                 = table(log(Y_rt(idx)),X_sal(idx),G_mou(idx),'VariableNames',{'RT','Saliency','Mouse'}); %Create table for mixed model
lme                 = fitlme(tbl,'RT~Modality+Saliency+(1|Mouse)'); %construct linear mixed effects model with fixed effect of cohort and random intercept for different mice
stats               = dataset2table(anova(lme,'DFMethod','Satterthwaite')); %Perform ANOVA on model and output as matrix

fprintf('Fixed effect of modality on RT (MST) (%d trials): (Linear Mixed Model ANOVA)\n',sum(idx))
fprintf('(F(%d,%2.0f) = %1.2f, p=%1.2e, ANOVA)\n',stats{2,3},stats{2,4},stats{2,2},stats{2,5})

tempfile = fullfile(params.savedir,'SourceData_Fig1h_Modality_RT.xlsx');
writetable(tbl,tempfile)

%Cohort on visual RT (UST vs MST):
X_coh       = NaN(nTrials,1);
X_coh(ismember(trialData.session_ID,sessionData.session_ID(ismember(sessionData.Experiment,params.Experiments{1}))))         = 1;
X_coh(ismember(trialData.session_ID,sessionData.session_ID(ismember(sessionData.Experiment,params.Experiments{2}))))         = 2;
X_coh(ismember(trialData.session_ID,sessionData.session_ID(ismember(sessionData.Experiment,params.Experiments{3}))))         = 3;

idx_vis             = strcmp(trialData.trialType,'X') & ismember(trialData.visualOriChangeNorm,2:5) & trialData.vecResponse==2 & idx_box & idx_lat;

idx                 = idx_vis & ismember(X_coh,[2 3]);

tbl                 = table(Y_rt(idx),X_coh(idx),X_sal_vis(idx),G_mou(idx),'VariableNames',{'RT','Cohort','Saliency','Mouse'}); %Create table for mixed model
lme                 = fitlme(tbl,'RT~Cohort+Saliency'); %construct linear mixed effects model with fixed effect of cohort and random intercept for different mice

stats               = dataset2table(anova(lme,'DFMethod','Satterthwaite')); %Perform ANOVA on model and output as matrix

fprintf('Fixed effect of cohort on visual RT (UST vs MST) (%d trials): (Linear Mixed Model ANOVA)\n',sum(idx))
fprintf('(F(%d,%2.0f) = %1.2f, p=%1.2e, ANOVA)\n',stats{2,3},stats{2,4},stats{2,2},stats{2,5})

tempfile = fullfile(params.savedir,'SourceData_Fig1h_Cohort_VisualRT.xlsx');
writetable(tbl,tempfile)

%% OLD STATS:
% %% Auditory vs visual RT in MST:
% p       = ranksum(datatotest(ismember(groups,5:8)),datatotest(ismember(groups,9:12)));
% fprintf('Wilcoxon rank sum test (%d visual, %d auditory trials): p=%1.2e\n',sum(ismember(groups,5:8)),sum(ismember(groups,9:12)),p)
% 
% %% Linear correlation RT and saliency for vision UST, vision MST and audition MST:
% [r,p] = corr(groups(ismember(groups,1:4)),datatotest(ismember(groups,1:4)));
% fprintf('RT - vis saliency UST: Pearson correlation (n=%d trials): r=%1.3f,p=%1.2e\n',sum(ismember(groups,1:4)),r,p)
% [r,p] = corr(groups(ismember(groups,5:8)),datatotest(ismember(groups,5:8)));
% fprintf('RT - vis saliency MST: Pearson correlation (n=%d trials): r=%1.3f,p=%1.2e\n',sum(ismember(groups,5:8)),r,p)
% [r,p] = corr(groups(ismember(groups,9:12)),datatotest(ismember(groups,9:12)));
% fprintf('RT - aud saliency MST: Pearson correlation (n=%d visual trials): r=%1.3f,p=%1.2e\n',sum(ismember(groups,9:12)),r,p)

% for i = 1:4 %For each saliency:
%     p       = ranksum(datatotest(groups==i),datatotest(groups==i+4)) * 12;
%     sigstar([i i+4],p) %use sigstar function to identify significance
%     fprintf('Wilcoxon rank sum test RT UST vs MST (Saliency level %d): %d UST trials, %d MST trial: p=%1.2e\n',i,sum(groups==i),sum(groups==i+4),p)
% end

%% Compute indices for hits of each of the different saliencies
params.trialcolors = {};
params.triallabels = {'Vsub' 'Vthr' 'Vsup' 'Vmax' 'Asub' 'Athr' 'Asup' 'Amax' };
splits              = {};
for iSplit = 2:5 %visual trials of different intensities:
    splits{end+1}               = strcmp(trialData.trialType,'X') & trialData.visualOriChangeNorm==iSplit & trialData.vecResponse==2;
    params.trialcolors{end+1}   = [0 0 iSplit/5];
end

for iSplit = 2:5%audio trials of different intensities:
    splits{end+1}               = strcmp(trialData.trialType,'Y') & trialData.audioFreqChangeNorm==iSplit & trialData.vecResponse==1; %#ok<*AGROW>
    params.trialcolors{end+1}   = [iSplit/5 0 0];
end
nSplits             = length(splits);

%% Compute dprime for each saliency level separately from fit:
medianRT_all    = NaN(nSessions,2);
dMat_all        = NaN(nSessions,8);

for iSes = 1:nSessions
    sesid               = sessionData.session_ID(iSes);
    [~,temptrialData]   = MOL_getTempPerSes(sesid,sessionData,trialData);%Get the sessionData for each session individually:
    
    xvalsau             = unique(abs(temptrialData.audioOctChange(ismember(temptrialData.audioFreqChangeNorm,2:5))));
    xvalsvis            = unique(abs(temptrialData.visualOriChange(ismember(temptrialData.visualOriChangeNorm,2:5))));
    
    dmax                = FullParam(iSes,1:2);
    n                   = FullParam(iSes,3:4);
    s50                 = FullParam(iSes,5:6);
    ce                  = FullParam(iSes,7:8);
    
    try
        % Compute d-prime for every value of visual and auditory:
        dMat_all(iSes,5:8) = dmax(1) * (xvalsau.^n(1))./(xvalsau.^n(1) + s50(1)^n(1)); %Auditory
        dMat_all(iSes,1:4) = dmax(2) * (xvalsvis.^n(2))./(xvalsvis.^n(2) + s50(2)^n(2)); %Visual
%         dMat_all(iSes,1:4) = dmax(1) * (xvalsau.^n(1))./(xvalsau.^n(1) + s50(1)^n(1)); %Auditory
%         dMat_all(iSes,5:8) = dmax(2) * (xvalsvis.^n(2))./(xvalsvis.^n(2) + s50(2)^n(2)); %Visual
    catch ME 
    end
    
    for iSplit = 1:length(splits)
        medianRT_all(iSes,iSplit)   = nanmedian(temptrialData.responseLatency(splits{iSplit}(strcmp(trialData.session_ID,sesid))));
    end
end
params.RTxlims = [180 800];

%% Figure showing correlation between reaction time and dprime (visual):
figure; hold all; set(gcf,'color','w','units','normalized','Position', [0.1 0.5 .3 .4]);
mrkrsize = 15;

%Visually trained mice:
idx_ses     = ismember(sessionData.Experiment,[params.Experiments(2) params.Experiments{3}]); hold all;

binres      = 100e3;
binedges    = (200e3:binres:850e3+binres)-binres/2;
binx        = 200e3:binres:850e3;
nBins       = length(binx);

handles = [];
splitidx = 1:4;
for iSplit = splitidx %nSplits
    datatoplot = NaN(1,nBins);
    
    Xdata       = medianRT_all(idx_ses,iSplit);
    Ydata       = dMat_all(idx_ses,iSplit);
%     idx         = ~isnan(Xdata) & ~isnan(Ydata); Xdata = Xdata(idx); Ydata = Ydata(idx);
    
    for iBin = 1:nBins
        idx = Xdata>binedges(iBin) & Xdata<binedges(iBin+1);
        if sum(idx)>2
            datatoplot(iBin) = nanmean(Ydata(idx));
        end
    end
    handles(end+1) = scatter(Xdata*1e-3,Ydata,mrkrsize,params.trialcolors{iSplit},'filled');
end

%Figure make up:
legend(handles,params.triallabels(splitidx),'Location','NorthEast'); legend boxoff
set(gca,'XTick',[200 500 800],'YTick',[1.75 3.5])
plot([0 1e6],[0 0],'k:','LineWidth',1)
xlim(params.RTxlims); ylim([0 3.5])
ylabel('Dprime')
xlabel('Median RT (ms)')

%Statistics:
G_mou       = cell(nSessions,1);
uMice       = unique(sessionData.mousename);
for iMouse = 1:length(uMice)
    G_mou(strcmp(sessionData.mousename,uMice{iMouse})) = uMice(iMouse);
end

X_RT                = reshape(medianRT_all(idx_ses,splitidx),numel(medianRT_all(idx_ses,splitidx)),1);
Y_Dprime            = reshape(dMat_all(idx_ses,splitidx),numel(dMat_all(idx_ses,splitidx)),1);
G_mou2              = repmat(G_mou(idx_ses),numel(splitidx),1);
X_sal               = reshape(repmat(1:4,sum(idx_ses),1),size(X_RT));

tbl                 = table(Y_Dprime,X_RT,X_sal,G_mou2,'VariableNames',{'Dprime','RT','Saliency','Mouse'}); %Create table for mixed model
lme                 = fitlme(tbl,'Dprime~RT+Saliency+(1|Mouse)'); %construct linear mixed effects model with fixed effect of cohort and random intercept for different mice
stats               = dataset2table(anova(lme,'DFMethod','Satterthwaite')); %Perform ANOVA on model and output as matrix

fprintf('Fixed effect of visual RT on dprime (%d sessions, 4 saliencies each): (Linear Mixed Model ANOVA)\n',sum(idx_ses))
fprintf('(F(%d,%2.0f) = %1.2f, p=%1.2e, ANOVA)\n',stats{2,3},stats{2,4},stats{2,2},stats{2,5})

tempfile = fullfile(params.savedir,'SourceData_FigS2g_VisDprime_RT.xlsx');
writetable(tbl,tempfile)

% [r,p]       = corr(Xdata',Ydata');
% fprintf('%d sessions, %d conditions\n',sum(idx_ses),numel(Xdata))
% if p<0.05
%     Fit = polyfit(Xdata',Ydata',1);
%     plot(100:1:1000,polyval(lme,100e3:1e3:1000e3),'k','LineWidth',2)
%     text(700,1.2,'***','Color','k','FontSize',50)
%     text(300,4-iSplit*0.25,sprintf('r=%1.2f,p=%1.3e',r,p),'Color','k','FontSize',15)
% end

%% Figure showing correlation between RT and dprime (auditory):
figure; hold all; set(gcf,'color','w','units','normalized','Position', [0.1 0.5 .3 .4]);
mrkrsize = 15;

%Multisensory trained mice (MST)
idx_ses     = ismember(sessionData.Experiment,params.Experiments{3}); hold all;

binres      = 100e3;
binedges    = (200e3:binres:850e3+binres)-binres/2;
binx        = 200e3:binres:850e3;
nBins       = length(binx);

handles = [];
splitidx = 5:8;
for iSplit = splitidx %nSplits
    datatoplot = NaN(1,nBins);
    
    Xdata       = medianRT_all(idx_ses,iSplit);
    Ydata       = dMat_all(idx_ses,iSplit);
    idx         = ~isnan(Xdata) & ~isnan(Ydata); Xdata = Xdata(idx); Ydata = Ydata(idx);
    
    for iBin = 1:nBins
        idx = Xdata>binedges(iBin) & Xdata<binedges(iBin+1);
        if sum(idx)>2
            datatoplot(iBin) = nanmean(Ydata(idx));
        end
    end
    handles(end+1) = scatter(Xdata*1e-3,Ydata,mrkrsize,params.trialcolors{iSplit},'filled');
end

%Figure make up:
legend(handles,params.triallabels(splitidx),'Location','NorthEast'); legend boxoff
set(gca,'XTick',[200 500 800],'YTick',[1.75 3.5])
plot([0 1e6],[0 0],'k:','LineWidth',1)
xlim(params.RTxlims); ylim([0 3.5])
ylabel('Dprime')
xlabel('Median RT (ms)')

X_RT                = reshape(medianRT_all(idx_ses,splitidx),numel(medianRT_all(idx_ses,splitidx)),1);
Y_Dprime            = reshape(dMat_all(idx_ses,splitidx),numel(dMat_all(idx_ses,splitidx)),1);
G_mou2              = repmat(G_mou(idx_ses),numel(splitidx),1);
X_sal               = reshape(repmat(1:4,sum(idx_ses),1),size(X_RT));

tbl                 = table(Y_Dprime,X_RT,X_sal,G_mou2,'VariableNames',{'Dprime','RT','Saliency','Mouse'}); %Create table for mixed model
lme                 = fitlme(tbl,'Dprime~RT+Saliency+(1|Mouse)'); %construct linear mixed effects model with fixed effect of cohort and random intercept for different mice
stats               = dataset2table(anova(lme,'DFMethod','Satterthwaite')); %Perform ANOVA on model and output as matrix

fprintf('Fixed effect of auditory RT on dprime (%d sessions, 4 saliencies each): (Linear Mixed Model ANOVA)\n',sum(idx_ses))
fprintf('(F(%d,%2.0f) = %1.2f, p=%1.2e, ANOVA)\n',stats{2,3},stats{2,4},stats{2,2},stats{2,5})

tempfile = fullfile(params.savedir,'SourceData_FigS2f_AudDprime_RT.xlsx');
writetable(tbl,tempfile)

% Xdata       = reshape(medianRT_all(:,splitidx),1,numel(medianRT_all(:,splitidx)));
% Ydata       = reshape(dMat_all(:,splitidx),1,numel(dMat_all(:,splitidx)));
% idx         = ~isnan(Xdata) & ~isnan(Ydata); Xdata = Xdata(idx); Ydata = Ydata(idx);
% 
% [r,p]       = corr(Xdata',Ydata');
% 
% [r,p]       = corr(Xdata',Ydata');
% fprintf('%d sessions, %d conditions\n',sum(idx_ses),numel(Xdata))
% if p<0.05
%     Fit = polyfit(Xdata',Ydata',1);
%     plot(100:1:1000,polyval(Fit,100e3:1e3:1000e3),'k','LineWidth',2)
%     text(700,1.2,'***','Color','k','FontSize',50)
%     text(300,4-iSplit*0.25,sprintf('r=%1.2f,p=%1.3e',r,p),'Color','k','FontSize',15)
% end

%% Cumulative distribution: figure not in manuscript

binres          = 25e3;

idx_lat         = ~(strcmp(trialData.trialType,'X') & trialData.responseLatency<200e3) | isnan(trialData.responseLatency);
idx_box         = ismember(trialData.session_ID,unique(sessionData.session_ID(sessionData.Box==9)));

iExp            = 2;
idx_exp         = ismember(trialData.session_ID,sessionData.session_ID(ismember(sessionData.Experiment,params.Experiments{iExp})));
allRT_UST       = trialData.responseLatency(strcmp(trialData.trialType,'X') & trialData.vecResponse==2 & idx_exp & idx_box & idx_lat);

iExp            = 3;
idx_exp         = ismember(trialData.session_ID,sessionData.session_ID(ismember(sessionData.Experiment,params.Experiments{iExp})));
allRT_MST       = trialData.responseLatency(strcmp(trialData.trialType,'X') & trialData.vecResponse==2 & idx_exp & idx_box & idx_lat);

%Convert to ms:

edges           = (0:binres:1e6) + binres/2; %#ok<NBRAK>
histmat         = [];
histmat(1,:)    = histcounts(allRT_UST,edges,'Normalization','cdf');
histmat(2,:)    = histcounts(allRT_MST,edges,'Normalization','cdf');

figure; set(gcf,'color','w','units','normalized','Position', [0.05 0.5 .35 .3]); hold all;
for iExp = 2:3
    plot(edges(1:end-1)+binres/2,squeeze(histmat(iExp-1,:)),'Color',params.colors_experiments{iExp})
end

%Figure makeup:
ylabel('Probability','FontSize',20)
xlabel('Reaction time (ms)','FontSize',20)
set(gca,'XTick',(0:200:1000)*1e3,'XTickLabel',0:200:1000)
legend(params.ExperimentLabels(2:3),'Location','SouthEast'); legend boxoff

%Statistics:
[p] = ranksum(allRT_UST,allRT_MST);
if p<1e-4
    text(200*1e3,0.8,'***','FontSize',75)
end

%%




%% MOL_Fig_Behavior_Supplementary figures:


%% Psychometric curves for octave and hertz separately:
params.auVersions = {'Hz' 'Oct'};

for iExp = 1:params.nExperiments
    
    for iAu = 1:2
        
        figure; set(gcf,'color','w','units','normalized','Position', [0.1 0.1 .7 .7]); hold all;
        expanimals                 = unique(sessionData.mousename(ismember(sessionData.Experiment,params.Experiments{iExp}) & ismember(sessionData.auChangeUnit,params.auVersions{iAu})));
        title(sprintf('%s, %s, (n=%d)',params.ExperimentLabels{iExp},params.auVersions{iAu},length(expanimals)))

        for iAnimal = 1:length(expanimals)
            
            [tempsessionData,~]         = MOL_getTempPerSes(sessionData.session_ID(strcmp(sessionData.mousename,expanimals{iAnimal})),sessionData,trialData);
            
            idx_ses                     = strcmp(sessionData.mousename,expanimals{iAnimal});
            meantheta_est               = median(FullParam(idx_ses,:),1);
            
            % Settings:
            switch tempsessionData.auChangeUnit{1}
                case 'Hz'
                    params.auprobepos       = 0.5;
                    params.auticks          = [20 100 4000];
                    params.auticklabels     = ['Probe' num2cell(params.auticks)];
                    params.auxaxislabel     = 'Delta frequency (Hz)';
                    params.auystatslabel    = 'Auditory threshold (Hz)';
                case 'Oct'
                    params.auprobepos       = 0.001;
                    params.auticks          = [1/64 1/8 1/2];
                    params.auticklabels     = {'Probe' '1/256' '1/64' '1/8' '1/2'};
                    params.auxaxislabel     = 'Delta frequency (Oct)';
                    params.auystatslabel    = 'Auditory threshold (partial octave)';
            end
            
            params.visprobepos     = 0.5;
            params.visticks        = [5 15 90];
            params.vistickslabels  = ['Probe' num2cell(params.visticks)];
            params.visxaxislabel   = 'Delta orientation (Degrees)';
            params.visystatslabel  = 'Visual threshold (Degrees)';
            
            params.yticks          = [0 0.25 0.5 0.75 1];
            
            % Generate contingency table from fitted parameters:
            [xvals_fit_au,xvals_fit_vis,ctable_fit_mat] = MOL_Gen2ADC_PsyCurve(meantheta_est,params);
            
            %         xvals_fit_au = 10.^(linspace(log10(0.001),log10(1/2),1000)); %logarithmic spacing
            
            %Audio:
            subplot(1,2,1); hold all;
%             plot(xvals_fit_au,squeeze(ctable_fit_mat(1,1,:)),'-','Color',params.cmap(iAnimal,:),'LineWidth',3);
%             plot(xvals_fit_au,squeeze(ctable_fit_mat(1,2,:)),':','Color',params.cmap(iAnimal,:),'LineWidth',3);
            
            plot(xvals_fit_au,squeeze(ctable_fit_mat(1,1,:)),'-','Color',[1 0 0],'LineWidth',3);
            plot(xvals_fit_au,squeeze(ctable_fit_mat(1,2,:)),':','Color',[0.2 0.2 1],'LineWidth',3);
            
            %Visual:
            subplot(1,2,2); hold all;
            plot(xvals_fit_vis,squeeze(ctable_fit_mat(2,2,:)),'-','Color',[0 0 1],'LineWidth',3);
            plot(xvals_fit_vis,squeeze(ctable_fit_mat(2,1,:)),':','Color',[1 0.2 0.2],'LineWidth',3);
        end
        
        MOL_Psy2Sided_FigMakeup(params)
        
        filename = sprintf('PsyCurves_%s_%s.eps',params.ExperimentLabels{iExp},params.auVersions{iAu});
        export_fig(fullfile(params.savedir,filename),gcf)
    end
end

%% Figures: dprimes and thresholds:

params.stat_labels      = {'DAudio' 'DVisual' 'AuditoryThreshold' 'VisualThreshold'};
params.stat_indices     = [1 2 5 6];
params.stat_upper       = [5 5 200 15; 5 5 1/8 15];
params.stat_lower       = [-0.2 -0.2 0 0; -0.2 -0.2 0 0];

params.cmap = hsv(11);
params.cmap = params.cmap([1 2 5 6 7 8 9 10 11],:); %remove yellow colors
%make figures of dprimes and thresholds:

for iExp = 1:params.nExperiments
    for iAu = 1:2
        figure; set(gcf,'color','w','units','normalized','Position', [0.1 0.3 .85 .28]); hold all;
        expanimals                  = unique(sessionData.mousename(ismember(sessionData.Experiment,params.Experiments{iExp}) & ismember(sessionData.auChangeUnit,params.auVersions{iAu})));
        for iStat = 1:4
            subplot(1,5,iStat); hold all;
            ST                          = params.stat_indices(iStat);
%             title(sprintf('%s, %s, %s, (n=%d)',params.stat_labels{iStat},params.ExperimentLabels{iExp},params.auVersions{iAu},length(expanimals)),'FontSize',15)
            title(sprintf('%s',params.stat_labels{iStat}),'FontSize',15)
            
            for iAnimal = 1:length(expanimals)
                idx_ses                     = strcmp(sessionData.mousename,expanimals{iAnimal});
                
                scatter(1+params.offset*iAnimal + randn(sum(idx_ses),1)*0.02,FullParam(idx_ses,ST),50,params.cmap(iAnimal,:),'filled')
            end
            
            idx_ses                     = ismember(sessionData.mousename,expanimals);

            boxplot(FullParam(idx_ses,ST),'plotstyle','compact','colors','k','medianstyle','target','whisker',1e6) %'k.','MarkerSize',40,'LineWidth',5)
            
            %Figure makeup:
%             ylabel(params.stat_labels{iStat})
            set(gca, 'XTick', 1,'XTickLabels', 'All','XTickLabelRotation',45)
%             ylim([0 params.stat_lims(iAu,iStat)])
            xlim([0.5 10.5])
            if iStat <3
                templim = get(gca,'Ylim');
                ylim([params.stat_lower(iAu,iStat) max([ceil(templim(2)) 2])])
            else
                ylim([params.stat_lower(iAu,iStat) params.stat_upper(iAu,iStat)])
            end
            
%             if iStat==3 && iAu==1
%                 set(gca, 'YScale', 'log')
%                 set(gca,'YTick',[10 100 500])
%                 ylim([1 params.stat_lims(iAu,iStat)])
%             end
%              if iStat==3 && iAu==2
%                 set(gca, 'YScale', 'log')
%                 set(gca,'YTick',[1/64 1/8])
%                 ylim([1/128 params.stat_lims(iAu,iStat)])
%             end

        end
        subplot(1,5,5); hold all;
        for iAnimal = 1:length(expanimals)
            scatter(1,2,50,params.cmap(iAnimal,:),'filled')
        end
        legend(expanimals,'FontSize',10)
        tightfig();
        
        filename = sprintf('%s%s.eps',params.ExperimentLabels{iExp},params.auVersions{iAu});
        export_fig(fullfile(params.savedir,filename),gcf)
    end
end

% distFig()

%% 


%Trim FullParam to number of animals:
FullParam = FullParam(:,1:length(mouseids),:);
%Compute some values:
FullParamReshape = reshape(FullParam,size(FullParam,1),size(FullParam,2)*size(FullParam,3));

%%Checks and balances:
Param_animalmean    = nanmean(FullParam,3);
%     Param_animalsem     = nanstd(FullParam,[],3) ./ sqrt(sum(~isnan(FullParam(1,:)))); %not necessary?
Param_fullmean      = nanmean(FullParamReshape,2);
Param_fullsem       = nanstd(FullParamReshape,[],2) / sqrt(sum(~isnan(FullParamReshape(1,:))));

% Dprime figure
figure; set(gcf,'color','w','units','normalized','Position', [0.1 0.2 .6 .5]);
subplot(1,2,1);
plot(1:length(mouseids),squeeze(FullParam(1,:,:)),'r.','MarkerSize',20); hold all;
plot(1:length(mouseids),Param_animalmean(1,:),'k.','MarkerSize',45); hold all;
errorbar(length(mouseids)+1,Param_fullmean(1),Param_fullsem(1),'-or','MarkerSize',10,'MarkerEdgeColor','red','MarkerFaceColor','red')

ylabel('d-Prime Auditory')
set(gca, 'XTick', 1:length(mouseids),'XTickLabels', mouseids,'XTickLabelRotation',45)
ylim([0 4])
xlim([0.5 length(mouseids)+1.5])

subplot(1,2,2);
plot(1:length(mouseids),squeeze(FullParam(2,:,:)),'b.','MarkerSize',20); hold all;
plot(1:length(mouseids),Param_animalmean(2,:),'k.','MarkerSize',45); hold all;
errorbar(length(mouseids)+1,Param_fullmean(2),Param_fullsem(2),'-ob','MarkerSize',10,'MarkerEdgeColor','blue','MarkerFaceColor','blue')

ylabel('d-Prime Visual')
set(gca, 'XTick', 1:length(mouseids),'XTickLabels', mouseids,'XTickLabelRotation',45)
ylim([0 4])
xlim([0.5 length(mouseids)+1.5])

% Threshold figure:
figure; set(gcf,'color','w','units','normalized','Position', [0.1 0.2 .6 .5]);
subplot(1,2,1);
plot(1:length(mouseids),squeeze(FullParam(5,:,:)),'r.','MarkerSize',20); hold all;
plot(1:length(mouseids),Param_animalmean(5,:),'k.','MarkerSize',45); hold all;
errorbar(length(mouseids)+1,Param_fullmean(5),Param_fullsem(5),'-or','MarkerSize',10,'MarkerEdgeColor','red','MarkerFaceColor','red')

ylabel(Par.auystatslabel)
set(gca, 'XTick', 1:length(mouseids),'XTickLabels', mouseids,'XTickLabelRotation',45)
switch sessionData.auChangeUnit{1}
    case 'Hz'
        set(gca, 'YScale', 'log')
        ylim([5 6000])
        set(gca,'YTick',[10 100 1000],'YTickLabels',[10 100 1000])
    case 'Oct'
        ylim([0 0.1])
end
xlim([0.5 length(mouseids)+1.5])

subplot(1,2,2);
plot(1:length(mouseids),squeeze(FullParam(6,:,:)),'b.','MarkerSize',20); hold all;
plot(1:length(mouseids),Param_animalmean(6,:),'k.','MarkerSize',45); hold all;
errorbar(length(mouseids)+1,Param_fullmean(6),Param_fullsem(6),'-ob','MarkerSize',10,'MarkerEdgeColor','blue','MarkerFaceColor','blue')

ylabel(Par.visystatslabel)
set(gca, 'XTick', 1:length(mouseids),'XTickLabels', mouseids,'XTickLabelRotation',45)
ylim([0 20])
xlim([0.5 length(mouseids)+1.5])
