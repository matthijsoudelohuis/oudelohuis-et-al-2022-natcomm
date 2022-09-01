%% Optogenetic silencing of V1 during audiovisual change detection
% This script analyzes the behavioral performance during V1 optogenetic silencing
% with early or late onset (0 versus 200 ms).
% Reproduces figure 5efgh
% Oude Lohuis et al. 2022 Nat Comms

%% load the data:
load('E:\Matlab\MOL_Analysis\AOudeLohuisetal_2022_NatComms\5Opto\Data5_3.mat','sessionData','trialData')
fprintf('Dataset: %d sessions, %d trials\n',length(sessionData.session_ID),length(trialData.session_ID));

%% General settings:
params.minPhotostimpower   = 2;    %Only include session in which power was >2mW total
params.sumTrialsCondition  = 30;   %Minimum number of summed trials across probe, Thr and Max within photostimcondition to compute dprime of

params                      = MOL_getColors_CHDET(params);

params.lines_visual_opto   = {'-' ':' ':'};
params.lines_audio_opto    = {'-' ':' ':'};

params.mrkrsize             = 50;

params.posthoctest          = 'bonferroni'; %Posthoc correction after Kruskal Wallis non parametric anova

%% Initialize structure for saving output fit parameters:
nSessions               = length(sessionData.session_ID);
dVis                    = NaN(nSessions,2,3); %Init matrix for storing all visual dprime data
dAud                    = NaN(nSessions,2,3); %Init matrix for storing all audio dprime data
TotalResp               = NaN(3,3,3,3,nSessions); %Init matrix for storing all hitrate data

cVis                    = NaN(nSessions,2,3); %Init matrix for storing all visual criterion data
cAud                    = NaN(nSessions,2,3); %Init matrix for storing all audio criterion data

%% Loop over sessions and fit 2ADC model to each session (separately for control, early and late silencing trials)
for iSes = 1:nSessions
    fprintf('\n\n Fitting session %d/%d\n\n',iSes,nSessions);
    sesid                                                       = sessionData.session_ID(iSes);
    [tempsessionData,temptrialData]                             = MOL_getTempPerSes(sesid,sessionData,trialData);%Get the sessionData for each session individually:
    
    trialFields = fieldnames(temptrialData);
    for iF = 1:length(trialFields)
        if any(strfind(trialFields{iF},'audio')) && iscell(temptrialData.(trialFields{iF}))
            temptrialData.(trialFields{iF}) = cell2vec(temptrialData.(trialFields{iF}))';
        end
    end
    
    [visconditions,auconditions,FullRespMat,FullnTrialsMat]     = MOL_Psy_GetTrialMatOpto(tempsessionData,temptrialData);
    %output is Au X Vis X response matrix (dimension 3 (response): layer 1 is fraction auditory, 2 visual, 3 no response)
    
    %Correction: If psychometric protocol, take intermediate values to compare with only 2 levels present
    if numel(visconditions) > 2
        FullRespMat = FullRespMat(:,[1 3 end],:,:);
        FullnTrialsMat = FullnTrialsMat(:,[1 3 end],:);
    end
    if numel(auconditions) > 2
        FullRespMat = FullRespMat([1 3 end],:,:,:);
        FullnTrialsMat = FullnTrialsMat([1 3 end],:,:);
    end
    
    %Compute d-prime for each condition:
    TotalResp(:,:,:,:,iSes)     = FullRespMat;
    FullRespMat_trialn          = FullRespMat .* repmat(FullnTrialsMat,1,1,1,3);
    
    for iTrial = 1:2
        for iOpto = 1:3
            outcome         = NaN(3,3);
            outcome(1,:)    = squeeze(FullRespMat_trialn(1,1+iTrial,iOpto,:)); %visual trials
            outcome(2,:)    = squeeze(FullRespMat_trialn(1+iTrial,1,iOpto,:)); %audio trials
            outcome(3,:)    = squeeze(FullRespMat_trialn(1,1,iOpto,:)); %probe trials
            outcome         = outcome(:,[2 1 3]); %Swap response coding, visual first.
            
            if ismember(tempsessionData.Experiment,{'VisOnlyTwolevels'}) %set NAN to 0 for auditory resposnes in UST mice:
                outcome(isnan(outcome)) = 0;
            end
            
            if nansum(outcome(:))>params.sumTrialsCondition
                [dVis(iSes,iTrial,iOpto),dAud(iSes,iTrial,iOpto),cVis(iSes,iTrial,iOpto),cAud(iSes,iTrial,iOpto)] = ...
                    MOL_Fit_2ADC_Full_Session(tempsessionData,trialData,0,outcome);
            end
        end
    end
end

fprintf('\n\n Finished fitting behavioral model.\n\n\n')

%%






%% Show average response rate figure for photoinhibition @V1 in MST animals:
idx_ses = strcmp(sessionData.PhotostimArea,'V1') & ...
    ismember(sessionData.Experiment,{'ChangeDetectionConflict'});

MOL_plotOptoBehavior_Rates(params,TotalResp(:,:,:,:,idx_ses))

%% Show average response rate figure for photoinhibition @S1 in MST animals:
idx_ses = strcmp(sessionData.PhotostimArea,'S1') & ...
    ismember(sessionData.Experiment,{'ChangeDetectionConflict'});

MOL_plotOptoBehavior_Rates(params,TotalResp(:,:,:,:,idx_ses))

%% Show average response rate figure for photoinhibition @V1 in UST animals:
idx_ses = strcmp(sessionData.PhotostimArea,'V1') & ...
    ismember(sessionData.Experiment,{'VisOnlyTwolevels'});

MOL_plotOptoBehavior_Rates(params,TotalResp(:,:,:,:,idx_ses))

%% Show average response rate figure for photoinhibition @V1 in UST animals:
idx_ses = strcmp(sessionData.PhotostimArea,'S1') & ...
    ismember(sessionData.Experiment,{'VisOnlyTwolevels'});

MOL_plotOptoBehavior_Rates(params,TotalResp(:,:,:,:,idx_ses))

%% Show average dprime figure for photoinhibition @V1:
idx_UST = strcmp(sessionData.PhotostimArea,'V1') & ...
    ismember(sessionData.Experiment,{'VisOnlyTwolevels'});

idx_MST = strcmp(sessionData.PhotostimArea,'V1') & ...
    ismember(sessionData.Experiment,{'ChangeDetectionConflict'});

Mice_UST = sessionData.mousename(idx_UST);
Mice_MST = sessionData.mousename(idx_MST);

MOL_plotOptoBehavior_Dprime(params,dVis(idx_UST,:,:),dAud(idx_UST,:,:),dVis(idx_MST,:,:),dAud(idx_MST,:,:),Mice_UST,Mice_MST)


%% Show average dprime figure for photoinhibition @S1:
idx_UST = strcmp(sessionData.PhotostimArea,'S1') & ...
    ismember(sessionData.Experiment,{'VisOnlyTwolevels'});

idx_MST = strcmp(sessionData.PhotostimArea,'S1') & ...
    ismember(sessionData.Experiment,{'ChangeDetectionConflict'});

Mice_UST = sessionData.mousename(idx_UST);
Mice_MST = sessionData.mousename(idx_MST);

MOL_plotOptoBehavior_Dprime(params,dVis(idx_UST,:,:),dAud(idx_UST,:,:),dVis(idx_MST,:,:),dAud(idx_MST,:,:),Mice_UST,Mice_MST)

%% Create Source Data file
G_Mou       = repmat(sessionData.mousename,1,2,3);
X_Sal       = cat(2,repmat({'Thr'},nSessions,1,3),repmat({'Max'},nSessions,1,3));
X_Opto      = cat(3,repmat({'Ctrl'},nSessions,2,1),repmat({'Early'},nSessions,2,1),repmat({'Late'},nSessions,2,1));
X_Cohort    = cell(nSessions,2,3);
X_Cohort(ismember(sessionData.Experiment,{'VisOnlyTwolevels'}),:,:) = {'UST'};
X_Cohort(ismember(sessionData.Experiment,{'ChangeDetectionConflict'}),:,:) = {'MST'};
X_OptoArea    = cell(nSessions,2,3);
X_OptoArea(strcmp(sessionData.PhotostimArea,'S1'),:,:) = {'S1'};
X_OptoArea(strcmp(sessionData.PhotostimArea,'V1'),:,:) = {'V1'};
tbl         = table(dVis(:),dAud(:),G_Mou(:),X_Sal(:),X_Opto(:),X_Cohort(:),X_OptoArea(:),...
    'VariableNames',{'Dprime_Vis','Dprime_Aud','Mouse','Saliency','Opto','Cohort','OptoArea'}); %Create table for mixed model

tempfile = fullfile('SourceData_Fig5_FigS6_Dprime.xlsx');
writetable(tbl,tempfile)


%% SUBSAMPLE ANALYSIS: small control analysis that checks if results for MST late inactivation effects hold if we subsample the #number of UST sessions from the MST data:
idx_UST = strcmp(sessionData.PhotostimArea,'V1') & ...
    ismember(sessionData.Experiment,{'VisOnlyTwolevels'});

idx_MST = strcmp(sessionData.PhotostimArea,'V1') & ...
    ismember(sessionData.Experiment,{'ChangeDetectionConflict'});

Mice_UST = sessionData.mousename(idx_UST);
Mice_MST = sessionData.mousename(idx_MST);

uMice_MST = unique(sessionData.mousename(idx_MST));
uMice_UST = unique(sessionData.mousename(idx_UST));

for iSub = 1:1000
    %Do stats on the combined dataset:
    temp        = find(idx_MST);
    idx_ses     = temp(randperm(length(temp),sum(idx_UST)));

    dVis_MST    = dVis(idx_ses,:,:);
    dMat        = squeeze(dVis_MST(:,1,[1 3]));
    Mice        = sessionData.mousename(idx_ses);
    
    %Statistical testing:    signed rank test with bonferroni correction:
    Y               = [dMat(:,1); dMat(:,2)];
    X_inac          = [ones(size(dMat,1),1); ones(size(dMat,1),1)*2];
    G_mou           = [Mice; Mice;];
    tbl             = table(Y,X_inac,G_mou,'VariableNames',{'Dprime','Inac','Mouse'}); %Create table for mixed model
    
    lme             = fitlme(tbl,'Dprime~Inac+(1|Mouse)'); %construct linear mixed effects model with fixed effect of temporal window and random intercept for different mice
    stats           = dataset2table(anova(lme,'DFMethod','Satterthwaite')); %Perform ANOVA on model and output as matrix
    pvals_thr(iSub)     = stats{2,5};

    dMat = squeeze(dVis_MST(:,2,[1 3]));
    
    %Statistical testing:    signed rank test with bonferroni correction:
    Y               = [dMat(:,1); dMat(:,2)];
    X_inac          = [ones(size(dMat,1),1); ones(size(dMat,1),1)*2];
    G_mou           = [Mice; Mice;];
    tbl             = table(Y,X_inac,G_mou,'VariableNames',{'Dprime','Inac','Mouse'}); %Create table for mixed model
    
    lme             = fitlme(tbl,'Dprime~Inac+(1|Mouse)'); %construct linear mixed effects model with fixed effect of temporal window and random intercept for different mice
    stats           = dataset2table(anova(lme,'DFMethod','Satterthwaite')); %Perform ANOVA on model and output as matrix
    pvals_max(iSub)     = stats{2,5};
    
end
sum(pvals_thr<0.05) / 1000
sum(pvals_max<0.05) / 1000


%% Compute median reaction time per session on visual trial types (Thr and Max change):
medianRT_all = NaN(nSessions,2);
for iSes = 1:nSessions
    sesid                           = sessionData.session_ID(iSes);
    [~,temptrialData]               = MOL_getTempPerSes(sesid,sessionData,trialData);%Get the sessionData for each session individually:
    medianRT_all(iSes,1)            = nanmedian(temptrialData.responseLatency(strcmp(temptrialData.trialType,'X') & temptrialData.visualOriChangeNorm==2 & temptrialData.hasphotostim~=1));
    medianRT_all(iSes,2)            = nanmedian(temptrialData.responseLatency(strcmp(temptrialData.trialType,'X') & temptrialData.visualOriChangeNorm==3 & temptrialData.hasphotostim~=1));
end

%% Compute reduction in dprime per session for Thr and Max changes and early and late inhibition:
delta_dVis                      = [];
delta_dVis(:,:,1)               = (dVis(:,:,1) - dVis(:,:,2)) ./ dVis(:,:,1) * 100; %d' early opto - d' control / d' control
delta_dVis(:,:,2)               = (dVis(:,:,1) - dVis(:,:,3)) ./ dVis(:,:,1) * 100; %d' late opto - d' control / d' control
delta_dVis(delta_dVis<-100)     = NaN; %Set to NaN, DECREASES in Dprime over 100% no sessions after updating fitting procedure
delta_dVis(delta_dVis>200)      = NaN; %Set to NaN, INCREASES in Dprime over 200% no sessions after updating fitting procedure

%% Figure that correlates the reduction in dprime early with late for MST mice:
% (not included in manuscript)
idx_ses = strcmp(sessionData.PhotostimArea,'V1') & ...
    ismember(sessionData.Experiment,{'ChangeDetectionConflict'});

figure; hold all; set(gcf,'color','w','units','normalized','Position', [0.1 0.4 .2 .3]);
X       = [delta_dVis(idx_ses,1,1); delta_dVis(idx_ses,2,1)];
Y       = [delta_dVis(idx_ses,1,2); delta_dVis(idx_ses,2,2)];
idx     = ~isnan(X) & ~isnan(Y);
X       = X(idx);
Y       = Y(idx);

%Figure make up:
lims = [-100 200];
scatter(X,Y,params.mrkrsize,params.colors_experiments{3},'filled');
xlabel('Early (Delta D'') ')
ylabel('Late (Delta D'') ')
xlim(lims); ylim(lims)
plot(lims,lims,'k--','LineWidth',1)

%test correlation:
[r,p]   = corr(X,Y);
text(-50,200,sprintf('corr: r=%1.2f,p=%1.3f',r,p),'FontSize',15)
if p<0.05
    text(-50,150,'*','FontSize',50)
end

%test significant difference:
[p]   = ranksum(X,Y);
text(-50,130,sprintf('ranksum: p=%1.3f',p),'FontSize',15)
if p<0.05
    text(-50,80,'*','FontSize',50)
end
    
%% Figure that correlates the reduction in dprime early with late for MST mice:
% (not included in manuscript)

idx_ses = strcmp(sessionData.PhotostimArea,'V1') & ...
    ismember(sessionData.Experiment,{'VisOnlyTwolevels'});

figure; hold all; set(gcf,'color','w','units','normalized','Position', [0.1 0.4 .2 .3]);
X       = [delta_dVis(idx_ses,1,1); delta_dVis(idx_ses,2,1)];
Y       = [delta_dVis(idx_ses,1,2); delta_dVis(idx_ses,2,2)];
idx     = ~isnan(X) & ~isnan(Y);
X       = X(idx);
Y       = Y(idx);

%Figure make up:
lims = [-100 200];
scatter(X,Y,params.mrkrsize,params.colors_experiments{2},'filled');
xlabel('Early (Delta D'') ')
ylabel('Late (Delta D'') ')
xlim(lims); ylim(lims)
plot(lims,lims,'k--','LineWidth',1)

%test correlation:
[r,p]   = corr(X,Y);
text(-50,200,sprintf('corr: r=%1.2f,p=%1.3f',r,p),'FontSize',15)
if p<0.05
    text(-50,150,'**','FontSize',50)
end

%test significant difference:
[p]   = ranksum(X,Y);
text(-50,130,sprintf('ranksum: p=%1.3f',p),'FontSize',15)
if p<0.05
    text(-50,80,'**','FontSize',50)
end

%% Figure: validation of relationship between dprime and reaction time in the subset of sessions with optical silencing (control trials only)
% (not included in manuscript)

figure; hold all; set(gcf,'color','w','units','normalized','Position', [0.1 0.5 .3 .4]);
mrkrsize = 80;
params.RTxlims = [250 800];

%Multisensory trained mice (MST)
idx_ses = ismember(sessionData.Experiment,{'ChangeDetectionConflict'});

Xdata       = medianRT_all(idx_ses,1);
Ydata       = dVis(idx_ses,1,1);
idx         = ~isnan(Xdata) & ~isnan(Ydata); Xdata = Xdata(idx); Ydata = Ydata(idx); 
scatter(Xdata*1e-3,Ydata,mrkrsize,params.colors_experiments{3},'filled'); %Thr change, control

Xdata       = medianRT_all(idx_ses,2);
Ydata       = dVis(idx_ses,2,1);
idx         = ~isnan(Xdata) & ~isnan(Ydata); Xdata = Xdata(idx); Ydata = Ydata(idx);
scatter(Xdata*1e-3,Ydata,mrkrsize,params.colors_experiments{3}.^2,'filled'); %Max change, control

Xdata       = medianRT_all(idx_ses,:); Xdata = Xdata(:);
Ydata       = dVis(idx_ses,:,1); Ydata  = Ydata(:);
idx         = ~isnan(Xdata) & ~isnan(Ydata); Xdata = Xdata(idx); Ydata = Ydata(idx);

[r,p]       = corr(Xdata,Ydata);
if p<0.05
    text(700,1.2,'***','Color','k','FontSize',50)
    text(660,1.55,sprintf('MST r=%1.2f,p=%1.3f',r,p),'Color','k','FontSize',15)
end

%Unisensory trained mice (UST)
idx_ses =    ismember(sessionData.Experiment,{'VisOnlyTwolevels'});
Xdata       = medianRT_all(idx_ses,1);
Ydata       = dVis(idx_ses,1,1);
idx         = ~isnan(Xdata) & ~isnan(Ydata); Xdata = Xdata(idx); Ydata = Ydata(idx); 
scatter(Xdata*1e-3,Ydata,mrkrsize,params.colors_experiments{2},'filled'); %Thr change, control

Xdata       = medianRT_all(idx_ses,2);
Ydata       = dVis(idx_ses,2,1);
idx         = ~isnan(Xdata) & ~isnan(Ydata); Xdata = Xdata(idx); Ydata = Ydata(idx);
scatter(Xdata*1e-3,Ydata,mrkrsize,params.colors_experiments{2}.^2,'filled'); %Max change, control

Xdata       = medianRT_all(idx_ses,:); Xdata = Xdata(:);
Ydata       = dVis(idx_ses,:,1); Ydata  = Ydata(:);
idx         = ~isnan(Xdata) & ~isnan(Ydata); Xdata = Xdata(idx); Ydata = Ydata(idx);

[r,p]       = corr(Xdata,Ydata);
if p<0.05
    text(700,0.7,'*','Color','k','FontSize',50)
    text(660,1,sprintf('UST r=%1.2f,p=%1.3f',r,p),'Color','k','FontSize',15)
end

legend({'MST Thr Change' 'MST Max Change' 'UST Thr Change' 'UST Max Change'}); legend boxoff
grid on;
plot([0 1e6],[0 0],'k:','LineWidth',1)
xlim(params.RTxlims); ylim([-0.5 3.5])
ylabel('Dprime')
xlabel('Median RT (ms)')

%% Compute reduction in dprime per session for Thr and Max changes and early and late inhibition:
delta_dVis          = [];
delta_dVis(:,:,1)   = (dVis(:,:,1) - dVis(:,:,2)) ./ dVis(:,:,1) * 100; %d' early opto - d' control / d' control
delta_dVis(:,:,2)   = (dVis(:,:,1) - dVis(:,:,3)) ./ dVis(:,:,1) * 100; %d' late opto - d' control / d' control
delta_dVis(delta_dVis<-100) = NaN; %Set to NaN, INCREASES in Dprime over 100% no sessions after updating fitting procedure
delta_dVis(delta_dVis>200) = NaN; %Set to NaN, DECREASES in Dprime over 200% no sessions after updating fitting procedure

%% For the coming figures, select only those sessions that have effective early silencing:
params.dprimereducthr   = 50; %In percentage decrease relative to control
% params.dprimereducthr   = 75; %In percentage decrease relative to control %uncomment to vary threshold
idx                     = delta_dVis(:,2,1)<params.dprimereducthr; %early silencing, Max change
delta_dVis(idx,:,:)     = NaN;
dVis_filter             = dVis;
dVis_filter(idx,:,:)    = NaN;

%% Figure that shows the dprime in control and either early or late silencing for each session connected by a line, as a function of median RT in that session
params.RTxlims = [250 600];

figure; hold all; set(gcf,'color','w','units','normalized','Position', [0.1 0.5 .64 .38]);
subplot(1,2,1); hold all;
title('Early Silencing')

mrkrshape = 'o';
mrkrsize = 40;

%Plot MST mice:
idx_ses = strcmp(sessionData.PhotostimArea,'V1') & ...
    ismember(sessionData.Experiment,{'ChangeDetectionConflict'});

scatter(medianRT_all(idx_ses,1)*1e-3,dVis_filter(idx_ses,1,1),mrkrsize,mrkrshape,'filled','MarkerFaceColor','w','MarkerEdgeColor',params.colors_experiments{3},'LineWidth',2); %Thr change, control - early
scatter(medianRT_all(idx_ses,2)*1e-3,dVis_filter(idx_ses,2,1),mrkrsize,mrkrshape,'filled','MarkerFaceColor',params.colors_experiments{3},'MarkerEdgeColor',params.colors_experiments{3},'LineWidth',2); %Max change, control - early
scatter(medianRT_all(idx_ses,1)*1e-3,dVis_filter(idx_ses,1,2),mrkrsize,'k','filled','MarkerEdgeColor','k','LineWidth',2); %Thr change, control - early
scatter(medianRT_all(idx_ses,2)*1e-3,dVis_filter(idx_ses,2,2),mrkrsize,'k','filled','MarkerEdgeColor','k','LineWidth',2); %Max change, control - early
plot([medianRT_all(idx_ses,1)*1e-3 medianRT_all(idx_ses,1)*1e-3]',[dVis_filter(idx_ses,1,1) dVis_filter(idx_ses,1,2)]','k-','LineWidth',1)
plot([medianRT_all(idx_ses,2)*1e-3 medianRT_all(idx_ses,2)*1e-3]',[dVis_filter(idx_ses,2,1) dVis_filter(idx_ses,2,2)]','k-','LineWidth',1)

%Plot UST mice:
idx_ses = strcmp(sessionData.PhotostimArea,'V1') & ...
    ismember(sessionData.Experiment,{'VisOnlyTwolevels'});

scatter(medianRT_all(idx_ses,1)*1e-3,dVis_filter(idx_ses,1,1),mrkrsize,mrkrshape,'filled','MarkerFaceColor','w','MarkerEdgeColor',params.colors_experiments{2},'LineWidth',2); %Thr change, control - early
scatter(medianRT_all(idx_ses,2)*1e-3,dVis_filter(idx_ses,2,1),mrkrsize,mrkrshape,'filled','MarkerFaceColor',params.colors_experiments{2},'MarkerEdgeColor',params.colors_experiments{2},'LineWidth',2); %Max change, control - early
scatter(medianRT_all(idx_ses,1)*1e-3,dVis_filter(idx_ses,1,2),mrkrsize,'k','filled','MarkerEdgeColor','k','LineWidth',2); %Thr change, control - early
scatter(medianRT_all(idx_ses,2)*1e-3,dVis_filter(idx_ses,2,2),mrkrsize,'k','filled','MarkerEdgeColor','k','LineWidth',2); %Max change, control - early
plot([medianRT_all(idx_ses,1)*1e-3 medianRT_all(idx_ses,1)*1e-3]',[dVis_filter(idx_ses,1,1) dVis_filter(idx_ses,1,2)]','k-','LineWidth',1)
plot([medianRT_all(idx_ses,2)*1e-3 medianRT_all(idx_ses,2)*1e-3]',[dVis_filter(idx_ses,2,1) dVis_filter(idx_ses,2,2)]','k-','LineWidth',1)

%Figure make up:
plot(params.RTxlims,[0 0],'k--','LineWidth',1)
xlim(params.RTxlims); ylim([-1.3 3.2])
ylabel('Dprime')
xlabel('Median RT (ms)')

subplot(1,2,2); hold all;
title('Late Silencing')

%Plot MST mice:
idx_ses = strcmp(sessionData.PhotostimArea,'V1') & ...
    ismember(sessionData.Experiment,{'ChangeDetectionConflict'});

scatter(medianRT_all(idx_ses,1)*1e-3,dVis_filter(idx_ses,1,1),mrkrsize,mrkrshape,'filled','MarkerFaceColor','w','MarkerEdgeColor',params.colors_experiments{3},'LineWidth',2); %Thr change, control - early
scatter(medianRT_all(idx_ses,2)*1e-3,dVis_filter(idx_ses,2,1),mrkrsize,mrkrshape,'filled','MarkerFaceColor',params.colors_experiments{3},'MarkerEdgeColor',params.colors_experiments{3},'LineWidth',2); %Max change, control - early
scatter(medianRT_all(idx_ses,1)*1e-3,dVis_filter(idx_ses,1,3),mrkrsize,'k','filled','MarkerEdgeColor','k','LineWidth',2); %Thr change, control - early
scatter(medianRT_all(idx_ses,2)*1e-3,dVis_filter(idx_ses,2,3),mrkrsize,'k','filled','MarkerEdgeColor','k','LineWidth',2); %Max change, control - early
plot([medianRT_all(idx_ses,1)*1e-3 medianRT_all(idx_ses,1)*1e-3]',[dVis_filter(idx_ses,1,1) dVis_filter(idx_ses,1,3)]','k-','LineWidth',1)
plot([medianRT_all(idx_ses,2)*1e-3 medianRT_all(idx_ses,2)*1e-3]',[dVis_filter(idx_ses,2,1) dVis_filter(idx_ses,2,3)]','k-','LineWidth',1)

%Plot UST mice:
idx_ses = strcmp(sessionData.PhotostimArea,'V1') & ...
    ismember(sessionData.Experiment,{'VisOnlyTwolevels'});

scatter(medianRT_all(idx_ses,1)*1e-3,dVis_filter(idx_ses,1,1),mrkrsize,mrkrshape,'filled','MarkerFaceColor','w','MarkerEdgeColor',params.colors_experiments{2},'LineWidth',2); %Thr change, control - early
scatter(medianRT_all(idx_ses,2)*1e-3,dVis_filter(idx_ses,2,1),mrkrsize,mrkrshape,'filled','MarkerFaceColor',params.colors_experiments{2},'MarkerEdgeColor',params.colors_experiments{2},'LineWidth',2); %Max change, control - early
scatter(medianRT_all(idx_ses,1)*1e-3,dVis_filter(idx_ses,1,3),mrkrsize,'k','filled','MarkerEdgeColor','k','LineWidth',2); %Thr change, control - early
scatter(medianRT_all(idx_ses,2)*1e-3,dVis_filter(idx_ses,2,3),mrkrsize,'k','filled','MarkerEdgeColor','k','LineWidth',2); %Max change, control - early
plot([medianRT_all(idx_ses,1)*1e-3 medianRT_all(idx_ses,1)*1e-3]',[dVis_filter(idx_ses,1,1) dVis_filter(idx_ses,1,3)]','k-','LineWidth',1)
plot([medianRT_all(idx_ses,2)*1e-3 medianRT_all(idx_ses,2)*1e-3]',[dVis_filter(idx_ses,2,1) dVis_filter(idx_ses,2,3)]','k-','LineWidth',1)

%Figure make up:
plot(params.RTxlims,[0 0],'k--','LineWidth',1)
xlim(params.RTxlims); ylim([-1.3 3.2])
ylabel('Dprime')
xlabel('Median RT (ms)')

%% Figure of relationship between reaction time and the effect of silencing on behavior: 

figure; hold all; set(gcf,'color','w','units','normalized','Position', [0.1 0.1 .3 .4]);
title('Early Silencing')

%Plot MST mice:
idx_ses = strcmp(sessionData.PhotostimArea,'V1') & ...
    ismember(sessionData.Experiment,{'ChangeDetectionConflict'});
scatter(medianRT_all(idx_ses,1)*1e-3,delta_dVis(idx_ses,1,1),mrkrsize,mrkrshape,'filled','MarkerFaceColor','w','MarkerEdgeColor',params.colors_experiments{3},'LineWidth',2); %Thr change, control - early
scatter(medianRT_all(idx_ses,2)*1e-3,delta_dVis(idx_ses,2,1),mrkrsize,mrkrshape,'filled','MarkerFaceColor',params.colors_experiments{3},'MarkerEdgeColor',params.colors_experiments{3},'LineWidth',2); %Max change, control - early

%Plot UST mice:
idx_ses = strcmp(sessionData.PhotostimArea,'V1') & ...
    ismember(sessionData.Experiment,{'VisOnlyTwolevels'});
scatter(medianRT_all(idx_ses,1)*1e-3,delta_dVis(idx_ses,1,1),mrkrsize,mrkrshape,'filled','MarkerFaceColor','w','MarkerEdgeColor',params.colors_experiments{2},'LineWidth',2); %Thr change, control - early
scatter(medianRT_all(idx_ses,2)*1e-3,delta_dVis(idx_ses,2,1),mrkrsize,mrkrshape,'filled','MarkerFaceColor',params.colors_experiments{2},'MarkerEdgeColor',params.colors_experiments{2},'LineWidth',2); %Max change, control - early

%Do stats on the combined dataset:
idx_ses     = strcmp(sessionData.PhotostimArea,'V1');
%Fit linear mixed model:
idx_thr         = ~isnan(delta_dVis(:,1,1)) & idx_ses;
idx_max         = ~isnan(delta_dVis(:,2,1)) & idx_ses;
X               = [medianRT_all(idx_thr,1); medianRT_all(idx_max,2)];
Y               = [delta_dVis(idx_thr,1,1); delta_dVis(idx_max,2,1)];
G_mou           = [sessionData.mousename(idx_thr); sessionData.mousename(idx_max)];
tbl             = table(Y,X*1e-3,G_mou,'VariableNames',{'Dprimereduc','RT','Mouse'}); %Create table for mixed model
lme             = fitlme(tbl,'Dprimereduc~RT+(1|Mouse)'); %construct linear mixed effects model with fixed effect of reaction time and random intercept for different mice
stats           = dataset2table(anova(lme,'DFMethod','Satterthwaite')); %Perform ANOVA on model and output as matrix
fprintf('Late inactivation:')
fprintf('F(%d,%2.0f)=%1.2f, p=%1.2f; \n',stats{2,3},stats{2,4},stats{2,2},stats{2,5})
text(550,150,sprintf('p=%1.3f',stats{2,5}),'Color','k','FontSize',15)
xvals = linspace(0,800,100);
plot(xvals,xvals*lme.Coefficients{2,2} + lme.Coefficients{1,2},'k-','LineWidth',5);

fprintf('%d THR %d MAX conditions, \n',sum(idx_thr), sum(idx_max))
title('')

tempfile = fullfile('SourceData_Fig5i_Dprime_RT_Early.xlsx');
writetable(tbl,tempfile)

%Figure make up:
% legend({'MST Thr Change' 'MST Max Change' 'UST Thr Change' 'UST Max Change'},'FontSize',15,'Location','SouthEast'); legend boxoff
plot(params.RTxlims,[0 0],'k--','LineWidth',1)
plot(params.RTxlims,[100 100],'k--','LineWidth',1)
xlim(params.RTxlims); ylim([-100 200])
ylabel('Delta Dprime')
xlabel('Median RT (ms)')

%%%%%%%% Do the same for late silencing:
figure; hold all; set(gcf,'color','w','units','normalized','Position', [0.4 0.1 .3 .4]);
title('Late Silencing')

%Plot MST mice:
idx_ses = strcmp(sessionData.PhotostimArea,'V1') & ...
    ismember(sessionData.Experiment,{'ChangeDetectionConflict'});
scatter(medianRT_all(idx_ses,1)*1e-3,delta_dVis(idx_ses,1,2),mrkrsize,mrkrshape,'filled','MarkerFaceColor','w','MarkerEdgeColor',params.colors_experiments{3},'LineWidth',2); %Thr change, control - early
scatter(medianRT_all(idx_ses,2)*1e-3,delta_dVis(idx_ses,2,2),mrkrsize,mrkrshape,'filled','MarkerFaceColor',params.colors_experiments{3},'MarkerEdgeColor',params.colors_experiments{3},'LineWidth',2); %Max change, control - early

%Plot UST mice:
idx_ses = strcmp(sessionData.PhotostimArea,'V1') & ...
    ismember(sessionData.Experiment,{'VisOnlyTwolevels'});
scatter(medianRT_all(idx_ses,1)*1e-3,delta_dVis(idx_ses,1,2),mrkrsize,mrkrshape,'filled','MarkerFaceColor','w','MarkerEdgeColor',params.colors_experiments{2},'LineWidth',2);%Thr change, control - early
scatter(medianRT_all(idx_ses,2)*1e-3,delta_dVis(idx_ses,2,2),mrkrsize,mrkrshape,'filled','MarkerFaceColor',params.colors_experiments{2},'MarkerEdgeColor',params.colors_experiments{2},'LineWidth',2); %Max change, control - early

%Do stats on the combined dataset:
idx_ses     = strcmp(sessionData.PhotostimArea,'V1');

%Fit linear mixed model:
idx_thr         = ~isnan(delta_dVis(:,1,2)) & idx_ses;
idx_max         = ~isnan(delta_dVis(:,2,2)) & idx_ses;
X               = [medianRT_all(idx_thr,1); medianRT_all(idx_max,2)];
Y               = [delta_dVis(idx_thr,1,2); delta_dVis(idx_max,2,2)];
G_mou           = [sessionData.mousename(idx_thr); sessionData.mousename(idx_max)];

tbl             = table(Y,X*1e-3,G_mou,'VariableNames',{'Dprimereduc','RT','Mouse'}); %Create table for mixed model
lme             = fitlme(tbl,'Dprimereduc~RT+(1|Mouse)'); %construct linear mixed effects model with fixed effect of reaction time and random intercept for different mice
stats           = dataset2table(anova(lme,'DFMethod','Satterthwaite')); %Perform ANOVA on model and output as matrix
fprintf('Late inactivation:')
fprintf('F(%d,%2.0f)=%1.2f, p=%1.2f; \n',stats{2,3},stats{2,4},stats{2,2},stats{2,5})
text(550,150,sprintf('p=%1.3f',stats{2,5}),'Color','k','FontSize',15)
xvals = linspace(0,800,100);
plot(xvals,xvals*lme.Coefficients{2,2} + lme.Coefficients{1,2},'k-','LineWidth',5);

fprintf('%d THR %d MAX conditions, \n',sum(idx_thr), sum(idx_max))
title('')

tempfile = fullfile('SourceData_Fig5j_Dprime_RT_Late.xlsx');
writetable(tbl,tempfile)

%Figure make up:
legend({'MST Thr Change' 'MST Max Change' 'UST Thr Change' 'UST Max Change'},'FontSize',15,'Location','SouthEast'); legend boxoff
plot(params.RTxlims,[0 0],'k--','LineWidth',1)
plot(params.RTxlims,[100 100],'k--','LineWidth',1)
xlim(params.RTxlims); ylim([-100 200])
ylabel('Delta Dprime')
xlabel('Median RT (ms)')

%% Figure: example session behavioral hit rates:
% example_fast
idx = [29 37];
for i=idx
    MOL_plotOptoBehavior_Rates(params,TotalResp(:,:,:,:,i))
    text(1,0.8,sprintf('RT=%3.1fms',nanmean(medianRT_all(i,:))*1e-3),'FontSize',20)
    fprintf('%s\n',sessionData.Experiment{i})
end

%% example_slow
idx = [24 35];
for i=idx
    MOL_plotOptoBehavior_Rates(params,TotalResp(:,:,:,:,i))
    text(1,0.8,sprintf('RT=%3.1fms',nanmean(medianRT_all(i,:))*1e-3),'FontSize',20)
    fprintf('%s\n',sessionData.Experiment{i})
end

% %% Show average response rate figure for photoinhibition @V1 in MST animals for fast versus slow sessions separately
% temp                = nanmean(medianRT_all,2);
% params.cutoff(1)    = prctile(temp,50); %Get fastest percentile (fastest sessions) based on median RT of sessions with early silencing
% params.cutoff(2)    = prctile(temp,50); %Get slowest percentile
% 
% idx_ses = strcmp(sessionData.PhotostimArea,'V1') & ...
%     ismember(sessionData.Experiment,{'ChangeDetectionConflict'}) & temp>=params.cutoff(2);
% 
% MOL_plotOptoBehavior_Rates(params,TotalResp(:,:,:,:,idx_ses))
% title('Slow')
% idx_ses = strcmp(sessionData.PhotostimArea,'V1') & ...
%     ismember(sessionData.Experiment,{'ChangeDetectionConflict'}) & temp<=params.cutoff(1);
% 
% MOL_plotOptoBehavior_Rates(params,TotalResp(:,:,:,:,idx_ses))
% title('Fast')
% 
% %% Show average dprime figure for photoinhibition @V1:
% idx_UST = strcmp(sessionData.PhotostimArea,'V1') & ...
%     ismember(sessionData.Experiment,{'VisOnlyTwolevels'}) & temp>=params.cutoff(2);
% 
% idx_MST = strcmp(sessionData.PhotostimArea,'V1') & ...
%     ismember(sessionData.Experiment,{'ChangeDetectionConflict'}) & temp>=params.cutoff(2);
% 
% Mice_UST = sessionData.mousename(idx_UST);
% Mice_MST = sessionData.mousename(idx_MST);
% 
% MOL_plotOptoBehavior_Dprime(params,dVis(idx_UST,:,:),dAud(idx_UST,:,:),dVis(idx_MST,:,:),dAud(idx_MST,:,:),Mice_UST,Mice_MST)
% title('Slow')
% 
% idx_UST = strcmp(sessionData.PhotostimArea,'V1') & ...
%     ismember(sessionData.Experiment,{'VisOnlyTwolevels'}) & temp<=params.cutoff(1);
% 
% idx_MST = strcmp(sessionData.PhotostimArea,'V1') & ...
%     ismember(sessionData.Experiment,{'ChangeDetectionConflict'}) & temp<=params.cutoff(1);
% 
% Mice_UST = sessionData.mousename(idx_UST);
% Mice_MST = sessionData.mousename(idx_MST);
% 
% MOL_plotOptoBehavior_Dprime(params,dVis(idx_UST,:,:),dAud(idx_UST,:,:),dVis(idx_MST,:,:),dAud(idx_MST,:,:),Mice_UST,Mice_MST)
% title('Fast')

%% Figure: Reduction in dprime for fastest sessions versus reduction in dprime for slowest sessions:
%(not included in manuscript)
meantoplot = [];
errortoplot = [];

temp                = nanmean(medianRT_all,2);
params.cutoff(1)    = prctile(temp,50); %Get fastest percentile (fastest sessions) based on median RT of sessions with early silencing
params.cutoff(2)    = prctile(temp,50); %Get slowest percentile

figure; hold all; set(gcf,'color','w','units','normalized','Position', [0.1 0.1 .2 .3]);
subplot(1,2,1); hold all;
title('Early Silencing')
FastSes             = delta_dVis(cat(3,medianRT_all<=params.cutoff(1),false(nSessions,2)));
meantoplot(1)       = nanmean(FastSes);
errortoplot(1)      = nanstd(FastSes) / sqrt(numel(FastSes));
SlowSes             = delta_dVis(cat(3,medianRT_all>=params.cutoff(2),false(nSessions,2)));
meantoplot(2)       = nanmean(SlowSes);
errortoplot(2)      = nanstd(SlowSes) / sqrt(numel(SlowSes));
bar([1 2],meantoplot); errorbar([1 2],meantoplot,errortoplot,'k.','LineWidth',2)
ylabel('%Reduction d''')
ylim([0 110])
set(gca,'XTick',[1 2],'XTickLabel',{'RT (50% Fast)' 'RT (50% Slow)'},'XTickLabelRotation',45)

% Compute significance:
Y               = [FastSes; SlowSes];
X_inac          = [ones(size(FastSes,1),1); ones(size(SlowSes,1),1)*2];
G_mou           = [sessionData.mousename(medianRT_all(:,1)<=params.cutoff(1)); sessionData.mousename(medianRT_all(:,2)<=params.cutoff(1));...
                    sessionData.mousename(medianRT_all(:,1)>=params.cutoff(2)); sessionData.mousename(medianRT_all(:,2)>=params.cutoff(2))];
tbl             = table(Y,X_inac,G_mou,'VariableNames',{'Dprime','Inac','Mouse'}); %Create table for mixed model

lme             = fitlme(tbl,'Dprime~Inac+(1|Mouse)'); %construct linear mixed effects model with fixed effect of temporal window and random intercept for different mice
stats           = dataset2table(anova(lme,'DFMethod','Satterthwaite')); %Perform ANOVA on model and output as matrix
p               = stats{2,5};
    
% p = ranksum(FastSes,SlowSes);
fprintf('Significance dprime reduction fast vs slow sessions early silencing:\n  %d fast sessions vs %d slow sessions: %1.5f\n\n',sum(~isnan(FastSes)),sum(~isnan(SlowSes)),p)
fprintf('F(%d,%2.0f)=%1.2f, p=%1.2f; \n',stats{2,3},stats{2,4},stats{2,2},stats{2,5})
if p<0.05
    sigstar([1 2],p)
end

subplot(1,2,2); hold all;
title('Late Silencing')
FastSes             = delta_dVis(cat(3,false(nSessions,2),medianRT_all<=params.cutoff(1)));
meantoplot(1)       = nanmean(FastSes);
errortoplot(1)      = nanstd(FastSes) / sqrt(numel(FastSes));
SlowSes             = delta_dVis(cat(3,false(nSessions,2),medianRT_all>=params.cutoff(2)));
meantoplot(2)       = nanmean(SlowSes);
errortoplot(2)      = nanstd(SlowSes) / sqrt(numel(SlowSes));
bar([1 2],meantoplot); errorbar([1 2],meantoplot,errortoplot,'k.','LineWidth',2)
ylabel('%Reduction d''')
ylim([0 110])
set(gca,'XTick',[1 2],'XTickLabel',{'RT (50% Fast)' 'RT (50% Slow)'},'XTickLabelRotation',45)

% Compute significance:
Y               = [FastSes; SlowSes];
X_inac          = [ones(size(FastSes,1),1); ones(size(SlowSes,1),1)*2];
G_mou           = [sessionData.mousename(medianRT_all(:,1)<=params.cutoff(1)); sessionData.mousename(medianRT_all(:,2)<=params.cutoff(1));...
                    sessionData.mousename(medianRT_all(:,1)>=params.cutoff(2)); sessionData.mousename(medianRT_all(:,2)>=params.cutoff(2))];
tbl             = table(Y,X_inac,G_mou,'VariableNames',{'Dprime','Inac','Mouse'}); %Create table for mixed model

lme             = fitlme(tbl,'Dprime~Inac+(1|Mouse)'); %construct linear mixed effects model with fixed effect of temporal window and random intercept for different mice
stats           = dataset2table(anova(lme,'DFMethod','Satterthwaite')); %Perform ANOVA on model and output as matrix
p               = stats{2,5};

% Compute significance:
fprintf('Significance dprime reduction fast vs slow sessions late silencing:\n  %d fast sessions vs %d slow sessions: %1.5f\n\n',sum(~isnan(FastSes)),sum(~isnan(SlowSes)),p)
fprintf('F(%d,%2.0f)=%1.2f, p=%1.2f; \n',stats{2,3},stats{2,4},stats{2,2},stats{2,5})
if p<0.05
    sigstar([1 2],p)
end

%% Figure: Extension of RT by late silencing:
% (in supplementary figures:)
trial_idx                   = false(6,length(trialData.session_ID));

trial_idx(1,:)              = strcmp(trialData.trialType,'X') & trialData.hasphotostim~=1 & trialData.vecResponse==2 & trialData.visualOriChangeNorm==3;
trial_idx(2,:)              = strcmp(trialData.trialType,'X') & trialData.hasphotostim==1 & trialData.vecResponse==2 & trialData.visualOriChangeNorm==3 & trialData.PostChangeOptoStart==0;
trial_idx(3,:)              = strcmp(trialData.trialType,'X') & trialData.hasphotostim==1 & trialData.vecResponse==2 & trialData.visualOriChangeNorm==3 & trialData.PostChangeOptoStart==0.2;

trial_idx(4,:)              = strcmp(trialData.trialType,'X') & trialData.hasphotostim~=1 & trialData.vecResponse==2 & trialData.visualOriChangeNorm==2;
trial_idx(5,:)              = strcmp(trialData.trialType,'X') & trialData.hasphotostim==1 & trialData.vecResponse==2 & trialData.visualOriChangeNorm==2 & trialData.PostChangeOptoStart==0;
trial_idx(6,:)              = strcmp(trialData.trialType,'X') & trialData.hasphotostim==1 & trialData.vecResponse==2 & trialData.visualOriChangeNorm==2 & trialData.PostChangeOptoStart==0.2;

effect_idx                  = ismember(trialData.session_ID,sessionData.session_ID(delta_dVis(:,2,1)>50));

ses_idx                     = ismember(trialData.session_ID,sessionData.session_ID(strcmp(sessionData.PhotostimArea,'V1')));

datatoplot                  = NaN(6,length(trialData.session_ID)); %init storing variable
for i = 1:6
    temp = trialData.responseLatency(trial_idx(i,:) & effect_idx' & ses_idx');
    %     temp = trialData.responseLatency(trial_idx(i,:));
    datatoplot(i,1:numel(temp)) = temp;
end

meantoplot          = nanmedian(datatoplot,2);
errortoplot         = nanstd(datatoplot,[],2) / sqrt(sum(sum(~isnan(datatoplot)))/6);

figure; set(gcf,'color','w','units','normalized','Position', [0.1 0.1 .25 .35]); hold all;
idx = [1 3 4 6];
h = errorbar(1:length(idx),meantoplot(idx),errortoplot(idx),'k.','LineWidth',2);
h.MarkerSize = 50;
set(gca,'XTick',1:4,'XTickLabel',{'Max Change Ctrl' 'Max Change OptoLate' 'Thr Change Ctrl' 'Thr Change OptoLate'},'XTickLabelRotation',45,'FontSize',15)
set(gca,'YTick',[300 350 400 450 500]*1e3,'YTickLabel',[300 350 400 450 500])
ylabel('Reaction Time')

%Statistics: %Fit linear mixed model:
uMice               = unique(sessionData.mousename);
nTrials             = length(trialData.session_ID);
T_mou               = cell(nTrials,1);
for iMouse = 1:length(uMice)
    T_mou(ismember(trialData.session_ID,sessionData.session_ID(strcmp(sessionData.mousename,uMice{iMouse})))) = uMice(iMouse);
end

Type                    = []; %init storing variable
G_Mou                   = [];
Y_RT                    = [];

for i = 1:6
    idx         = trial_idx(i,:) & effect_idx' & ses_idx';
    Y_RT        = [Y_RT; trialData.responseLatency(idx)]; %#ok<*AGROW>
    Type        = [Type; ones(sum(idx),1)*i];
    G_Mou       = [G_Mou; T_mou(idx)];
end

idx             = ismember(Type,[1 3]);
tbl             = table(Y_RT(idx),Type(idx),G_Mou(idx),'VariableNames',{'RT','Inac','Mouse'}); %Create table for mixed model
lme             = fitglme(tbl,'RT~Inac+(1|Mouse)','Distribution','gamma','link','log'); %construct generalized linear mixed effects model with fixed effect of inactivation and random intercept for different mice
stats           = dataset2table(anova(lme,'DFMethod','Residual')); %Perform ANOVA on model and output as matrix
fprintf('Generalized Linear Mixed Model ANOVA test for difference of RT between:\n maximal change control hits (n=%d trials) and late silencing hits (n=%d trials):\n',sum(~isnan(datatoplot(1,:))),sum(~isnan(datatoplot(3,:))));
fprintf('F(%d,%2.0f)=%1.2f, p=%1.2f; \n',stats{2,3},stats{2,4},stats{2,2},stats{2,5})
sigstar([1 2],stats{2,5})

tempfile = fullfile('SourceData_Fig6l_RT_Inactivation_Early.xlsx');
writetable(tbl,tempfile)

idx             = ismember(Type,[4 6]);
tbl             = table(Y_RT(idx),Type(idx),G_Mou(idx),'VariableNames',{'RT','Inac','Mouse'}); %Create table for mixed model
lme             = fitglme(tbl,'RT~Inac+(1|Mouse)','Distribution','gamma','link','log'); %construct generalized linear mixed effects model with fixed effect of inactivation and random intercept for different mice
stats           = dataset2table(anova(lme,'DFMethod','Residual')); %Perform ANOVA on model and output as matrix
fprintf('Generalized Linear Mixed Model ANOVA test for difference of RT between:\n threshold change control hits (n=%d trials) and late silencing hits (n=%d trials):\n',sum(~isnan(datatoplot(4,:))),sum(~isnan(datatoplot(6,:))));
fprintf('F(%d,%2.0f)=%1.2f, p=%1.2f; \n',stats{2,3},stats{2,4},stats{2,2},stats{2,5})
sigstar([3 4],stats{2,5})

tempfile = fullfile('SourceData_Fig6l_RT_Inactivation_Late.xlsx');
writetable(tbl,tempfile)

%% Control figure: Figure of relationship between reaction time and the effect of silencing on criterion (not d'): 

%Compute reduction in dprime per session for Thr and Max changes and early and late inhibition:
delta_cVis          = [];
delta_cVis(:,:,1)   = (cVis(:,:,1) - cVis(:,:,2)) ./ cVis(:,:,1) * 100; %d' early opto - d' control / d' control
delta_cVis(:,:,2)   = (cVis(:,:,1) - cVis(:,:,3)) ./ cVis(:,:,1) * 100; %d' late opto - d' control / d' control

figure; hold all; set(gcf,'color','w','units','normalized','Position', [0.1 0.1 .3 .4]);
title('Early Silencing')

%Plot MST mice:
idx_ses = strcmp(sessionData.PhotostimArea,'V1') & ...
    ismember(sessionData.Experiment,{'ChangeDetectionConflict'});
scatter(medianRT_all(idx_ses,1)*1e-3,delta_cVis(idx_ses,1,1),mrkrsize,mrkrshape,'filled','MarkerFaceColor','w','MarkerEdgeColor',params.colors_experiments{3},'LineWidth',2); %Thr change, control - early
scatter(medianRT_all(idx_ses,2)*1e-3,delta_cVis(idx_ses,2,1),mrkrsize,mrkrshape,'filled','MarkerFaceColor',params.colors_experiments{3},'MarkerEdgeColor',params.colors_experiments{3},'LineWidth',2); %Max change, control - early

%Plot UST mice:
idx_ses = strcmp(sessionData.PhotostimArea,'V1') & ...
    ismember(sessionData.Experiment,{'VisOnlyTwolevels'});
scatter(medianRT_all(idx_ses,1)*1e-3,delta_cVis(idx_ses,1,1),mrkrsize,mrkrshape,'filled','MarkerFaceColor','w','MarkerEdgeColor',params.colors_experiments{2},'LineWidth',2);%Thr change, control - early
scatter(medianRT_all(idx_ses,2)*1e-3,delta_cVis(idx_ses,2,1),mrkrsize,mrkrshape,'filled','MarkerFaceColor',params.colors_experiments{2},'MarkerEdgeColor',params.colors_experiments{2},'LineWidth',2); %Max change, control - early

%Do stats on the combined dataset:
idx_ses     = strcmp(sessionData.PhotostimArea,'V1');

%Fit linear mixed model:
idx_thr         = ~isnan(delta_cVis(:,1,1)) & idx_ses;
idx_max         = ~isnan(delta_cVis(:,2,1)) & idx_ses;
X               = [medianRT_all(idx_thr,1); medianRT_all(idx_max,2)];
Y               = [delta_cVis(idx_thr,1,1); delta_cVis(idx_max,2,1)];
G_mou           = [sessionData.mousename(idx_thr); sessionData.mousename(idx_max)];

tbl             = table(Y,X*1e-3,G_mou,'VariableNames',{'Dprimereduc','RT','Mouse'}); %Create table for mixed model
lme             = fitlme(tbl,'Dprimereduc~RT+(1|Mouse)'); %construct linear mixed effects model with fixed effect of reaction time and random intercept for different mice
stats           = dataset2table(anova(lme,'DFMethod','Satterthwaite')); %Perform ANOVA on model and output as matrix
fprintf('Early inactivation:')
fprintf('F(%d,%2.0f)=%1.2f, p=%1.2f; \n',stats{2,3},stats{2,4},stats{2,2},stats{2,5})
text(550,150,sprintf('p=%1.3f',stats{2,5}),'Color','k','FontSize',15)
xvals = linspace(0,800,100);
plot(xvals,xvals*lme.Coefficients{2,2} + lme.Coefficients{1,2},'k-','LineWidth',5);

fprintf('%d THR %d MAX conditions, \n',sum(idx_thr), sum(idx_max))
title('')

%Figure make up:
legend({'MST Thr Change' 'MST Max Change' 'UST Thr Change' 'UST Max Change'},'FontSize',15,'Location','SouthEast'); legend boxoff
plot(params.RTxlims,[0 0],'k--','LineWidth',1)
plot(params.RTxlims,[100 100],'k--','LineWidth',1)
xlim(params.RTxlims); ylim([-200 200])
ylabel('Delta Criterion')
xlabel('Median RT (ms)')

%%%%%%%% Do the same for late silencing:
figure; hold all; set(gcf,'color','w','units','normalized','Position', [0.4 0.1 .3 .4]);
title('Late Silencing')

%Plot MST mice:
idx_ses = strcmp(sessionData.PhotostimArea,'V1') & ...
    ismember(sessionData.Experiment,{'ChangeDetectionConflict'});
scatter(medianRT_all(idx_ses,1)*1e-3,delta_cVis(idx_ses,1,2),mrkrsize,mrkrshape,'filled','MarkerFaceColor','w','MarkerEdgeColor',params.colors_experiments{3},'LineWidth',2); %Thr change, control - early
scatter(medianRT_all(idx_ses,2)*1e-3,delta_cVis(idx_ses,2,2),mrkrsize,mrkrshape,'filled','MarkerFaceColor',params.colors_experiments{3},'MarkerEdgeColor',params.colors_experiments{3},'LineWidth',2); %Max change, control - early

%Plot UST mice:
idx_ses = strcmp(sessionData.PhotostimArea,'V1') & ...
    ismember(sessionData.Experiment,{'VisOnlyTwolevels'});
scatter(medianRT_all(idx_ses,1)*1e-3,delta_cVis(idx_ses,1,2),mrkrsize,mrkrshape,'filled','MarkerFaceColor','w','MarkerEdgeColor',params.colors_experiments{2},'LineWidth',2);%Thr change, control - early
scatter(medianRT_all(idx_ses,2)*1e-3,delta_cVis(idx_ses,2,2),mrkrsize,mrkrshape,'filled','MarkerFaceColor',params.colors_experiments{2},'MarkerEdgeColor',params.colors_experiments{2},'LineWidth',2); %Max change, control - early

%Do stats on the combined dataset:
idx_ses     = strcmp(sessionData.PhotostimArea,'V1');

%Fit linear mixed model:
idx_thr         = ~isnan(delta_cVis(:,1,2)) & idx_ses;
idx_max         = ~isnan(delta_cVis(:,2,2)) & idx_ses;
X               = [medianRT_all(idx_thr,1); medianRT_all(idx_max,2)];
Y               = [delta_cVis(idx_thr,1,2); delta_cVis(idx_max,2,2)];
G_mou           = [sessionData.mousename(idx_thr); sessionData.mousename(idx_max)];

tbl             = table(Y,X*1e-3,G_mou,'VariableNames',{'Dprimereduc','RT','Mouse'}); %Create table for mixed model
lme             = fitlme(tbl,'Dprimereduc~RT+(1|Mouse)'); %construct linear mixed effects model with fixed effect of reaction time and random intercept for different mice
stats           = dataset2table(anova(lme,'DFMethod','Satterthwaite')); %Perform ANOVA on model and output as matrix
fprintf('Early inactivation:')
fprintf('F(%d,%2.0f)=%1.2f, p=%1.2f; \n',stats{2,3},stats{2,4},stats{2,2},stats{2,5})
text(550,150,sprintf('p=%1.3f',stats{2,5}),'Color','k','FontSize',15)
xvals = linspace(0,800,100);
plot(xvals,xvals*lme.Coefficients{2,2} + lme.Coefficients{1,2},'k-','LineWidth',5);

fprintf('%d THR %d MAX conditions, \n',sum(idx_thr), sum(idx_max))
title('')

tempfile = fullfile('SourceData_Fig6k_Criterion_RT_Late.xlsx');
writetable(tbl,tempfile)

%Figure make up:
legend({'MST Thr Change' 'MST Max Change' 'UST Thr Change' 'UST Max Change'},'FontSize',15,'Location','SouthEast'); legend boxoff
plot(params.RTxlims,[0 0],'k--','LineWidth',1)
plot(params.RTxlims,[100 100],'k--','LineWidth',1)
xlim(params.RTxlims); ylim([-200 200])
ylabel('Delta Criterion')
xlabel('Median RT (ms)')

%% Control figure: Figure of relationship between reaction time and FA rate to visual lick spout (not d'): 
figure; hold all; set(gcf,'color','w','units','normalized','Position', [0.4 0.1 .3 .4]);
title('FA Rate - Late Silencing')
delta_FA            = squeeze(TotalResp(1,1,3,2,:) - TotalResp(1,1,1,2,:)); %get change in visual hit rate (TotalResp is 5D is Aud X Vis X Opto X Resp X Session matrix)

%Plot MST mice:
idx_ses = strcmp(sessionData.PhotostimArea,'V1') & ...
    ismember(sessionData.Experiment,{'ChangeDetectionConflict'});
scatter(medianRT_all(idx_ses,2)*1e-3,delta_FA(idx_ses),mrkrsize,mrkrshape,'filled','MarkerFaceColor',params.colors_experiments{3},'MarkerEdgeColor',params.colors_experiments{3},'LineWidth',2); %Max change, control - early

%Plot UST mice:
idx_ses = strcmp(sessionData.PhotostimArea,'V1') & ...
    ismember(sessionData.Experiment,{'VisOnlyTwolevels'});
scatter(medianRT_all(idx_ses,2)*1e-3,delta_FA(idx_ses),mrkrsize,mrkrshape,'filled','MarkerFaceColor',params.colors_experiments{2},'MarkerEdgeColor',params.colors_experiments{2},'LineWidth',2); %Max change, control - early

%Do stats on the combined dataset:
idx_ses     = strcmp(sessionData.PhotostimArea,'V1');

%Fit linear mixed model:
idx             = ~isnan(delta_FA) & idx_ses;
X               = nanmean(medianRT_all,2);
Y               = delta_FA;
G_mou           = sessionData.mousename;

tbl             = table(Y(idx),X(idx)*1e-3,G_mou(idx),'VariableNames',{'Dprimereduc','RT','Mouse'}); %Create table for mixed model
lme             = fitlme(tbl,'Dprimereduc~RT+(1|Mouse)'); %construct linear mixed effects model with fixed effect of reaction time and random intercept for different mice
stats           = dataset2table(anova(lme,'DFMethod','Satterthwaite')); %Perform ANOVA on model and output as matrix
fprintf('Late inactivation:')
fprintf('F(%d,%2.0f)=%1.2f, p=%1.2f; \n',stats{2,3},stats{2,4},stats{2,2},stats{2,5})
text(550,150,sprintf('p=%1.3f',stats{2,5}),'Color','k','FontSize',15)
xvals = linspace(0,800,100);
plot(xvals,xvals*lme.Coefficients{2,2} + lme.Coefficients{1,2},'k-','LineWidth',5);

fprintf('%d THR %d MAX conditions, \n',sum(idx_thr), sum(idx_max))
title('')

%Figure make up:
legend({'MST Thr Change' 'MST Max Change' 'UST Thr Change' 'UST Max Change'},'FontSize',15,'Location','SouthEast'); legend boxoff
plot(params.RTxlims,[0 0],'k--','LineWidth',1)
plot(params.RTxlims,[100 100],'k--','LineWidth',1)
xlim(params.RTxlims); ylim([-0.65 0.65])
ylabel('Delta False alarm')
xlabel('Median RT (ms)')

%% Control analysis: Fit linear mixed model with visual dprime as covariate to account for the negative correlation between RT and visual Dprime
idx_thr         = ~isnan(delta_dVis(:,1,2)) & idx_ses;
idx_max         = ~isnan(delta_dVis(:,2,2)) & idx_ses;
X_dVis          = [dVis(idx_thr,1,1); dVis(idx_max,2,1)];
X_RT            = [medianRT_all(idx_thr,1); medianRT_all(idx_max,2)];
Y               = [delta_dVis(idx_thr,1,2); delta_dVis(idx_max,2,2)];
G_mou           = [sessionData.mousename(idx_thr); sessionData.mousename(idx_max)];

tbl             = table(Y,X_RT,X_dVis,G_mou,'VariableNames',{'Dprimereduc','RT','DVis','Mouse'}); %Create table for mixed model
lme             = fitlme(tbl,'Dprimereduc~RT+DVis+(1|Mouse)'); %construct linear mixed effects model with fixed effect of reaction time and random intercept for different mice
stats           = dataset2table(anova(lme,'DFMethod','Satterthwaite')); %Perform ANOVA on model and output as matrix
fprintf('Multivariate Linear Mixed Model (Late inactivation:)\n')
fprintf('RT: F(%d,%2.0f)=%1.2f, p=%1.3f; \n',stats{2,3},stats{2,4},stats{2,2},stats{2,5})
fprintf('DVis: F(%d,%2.0f)=%1.2f, p=%1.3f; \n',stats{3,3},stats{3,4},stats{3,2},stats{3,5})
fprintf('%d THR %d MAX conditions, \n',sum(idx_thr), sum(idx_max))


%DEPRECATED
% % %% Figure: Reduction in hit rate for 25% fastest sessions versus slowest 25% of sessions:
% % % % %(not included in manuscript)
% % meantoplot = [];
% % errortoplot = [];
% % 
% % temp                = nanmean(medianRT_all,2);
% % params.cutoff(1)    = prctile(temp,33); %Get 25th percentile (fastest sessions) based on median RT of sessions with early silencing
% % params.cutoff(2)    = prctile(temp,66); %Get 75th percentile
% % 
% % figure; hold all; set(gcf,'color','w','units','normalized','Position', [0.1 0.1 .2 .3]);
% % subplot(1,2,1); hold all;
% % title('Early Silencing')
% % 
% % delta_dHit          = [];
% % delta_dHit(:,1)     = TotalResp(1,2,1,2,:) - TotalResp(1,2,2,2,:); %v hit rate thr control - early
% % delta_dHit(:,2)     = TotalResp(1,3,1,2,:) - TotalResp(1,3,2,2,:); %v hit rate max control - early
% % 
% % % delta_dHit(:,1) = (TotalResp(1,2,1,2,:) - TotalResp(1,2,2,2,:)) ./ TotalResp(1,2,1,2,:); %v hit rate thr control - early
% % % delta_dHit(:,2) = (TotalResp(1,3,1,2,:) - TotalResp(1,3,2,2,:)) ./ TotalResp(1,3,1,2,:); %v hit rate max control - early
% % %     output is Au X Vis X response matrix (dimension 3 (response): layer 1 is fraction auditory, 2 visual, 3 no response)
% % 
% % % delta_dHit(isinf(delta_dHit))=NaN;
% % delta_dHit          = nanmean(delta_dHit,2)*1.5;
% % delta_dHit          = delta_dHit*100;
% % 
% % FastSes             = delta_dHit(temp<=params.cutoff(1));
% % meantoplot(1)       = nanmean(FastSes);
% % errortoplot(1)      = nanstd(FastSes) / sqrt(numel(FastSes));
% % SlowSes             = delta_dHit(temp>=params.cutoff(2));
% % meantoplot(2)       = nanmean(SlowSes);
% % errortoplot(2)      = nanstd(SlowSes) / sqrt(numel(SlowSes));
% % bar([1 2],meantoplot); errorbar([1 2],meantoplot,errortoplot,'k.','LineWidth',2)
% % ylabel('Reduction Hit rate (%)')
% % ylim([-5 30])
% % set(gca,'XTick',[1 2],'XTickLabel',{'RT (50% Fast)' 'RT (50% Slow)'},'XTickLabelRotation',45,'YTick',[0 10 20 30 40])
% % 
% % % Compute significance:
% % Y               = [FastSes; SlowSes];
% % X_inac          = [ones(size(FastSes,1),1); ones(size(SlowSes,1),1)*2];
% % G_mou           = [sessionData.mousename(temp<=params.cutoff(1)); sessionData.mousename(temp>=params.cutoff(2))];
% % tbl             = table(Y,X_inac,G_mou,'VariableNames',{'Dprime','Inac','Mouse'}); %Create table for mixed model
% % 
% % lme             = fitlme(tbl,'Dprime~Inac+(1|Mouse)'); %construct linear mixed effects model with fixed effect of temporal window and random intercept for different mice
% % stats           = dataset2table(anova(lme,'DFMethod','Satterthwaite')); %Perform ANOVA on model and output as matrix
% % p               = stats{2,5};
% % 
% % % p = ranksum(FastSes,SlowSes);
% % fprintf('Significance hit rate reduction difference (early silencing):\n  %d fast sessions vs %d slow sessions: %1.5f\n\n',sum(~isnan(FastSes)),sum(~isnan(SlowSes)),p)
% % fprintf('F(%d,%2.0f)=%1.2f, p=%1.2f; \n',stats{2,3},stats{2,4},stats{2,2},stats{2,5})
% % if p<0.05
% %     sigstar([1 2],p)
% % end
% % 
% % delta_dHit          = [];
% % delta_dHit(:,1)     = TotalResp(1,2,1,2,:) - TotalResp(1,2,3,2,:); %v hit rate thr control - late
% % delta_dHit(:,2)     = TotalResp(1,3,1,2,:) - TotalResp(1,3,3,2,:); %v hit rate max control - late
% % 
% % delta_dHit          = nanmean(delta_dHit,2)*1.5;
% % delta_dHit          = delta_dHit*100;
% % 
% % subplot(1,2,2); hold all;
% % title('Late Silencing')
% % FastSes             = delta_dHit(temp<=params.cutoff(1));
% % meantoplot(1)       = nanmean(FastSes);
% % errortoplot(1)      = nanstd(FastSes) / sqrt(numel(FastSes));
% % SlowSes             = delta_dHit(temp>=params.cutoff(2));
% % meantoplot(2)       = nanmean(SlowSes);
% % errortoplot(2)      = nanstd(SlowSes) / sqrt(numel(SlowSes));
% % bar([1 2],meantoplot); errorbar([1 2],meantoplot,errortoplot,'k.','LineWidth',2)
% % % ylabel('Reduction Hit rate (%)')
% % ylim([-2.5 15])
% % set(gca,'XTick',[1 2],'XTickLabel',{'RT (50% Fast)' 'RT (50% Slow)'},'XTickLabelRotation',45,'YTick',[0 5 10 15 20],'YTickLabels',[0 10 20 30 40])
% % 
% % % Compute significance:
% % Y               = [FastSes; SlowSes];
% % X_inac          = [ones(size(FastSes,1),1); ones(size(SlowSes,1),1)*2];
% % G_mou           = [sessionData.mousename(temp<=params.cutoff(1)); sessionData.mousename(temp>=params.cutoff(2))];
% % tbl             = table(Y,X_inac,G_mou,'VariableNames',{'Dprime','Inac','Mouse'}); %Create table for mixed model
% % 
% % lme             = fitlme(tbl,'Dprime~Inac+(1|Mouse)'); %construct linear mixed effects model with fixed effect of temporal window and random intercept for different mice
% % stats           = dataset2table(anova(lme,'DFMethod','Satterthwaite')); %Perform ANOVA on model and output as matrix
% % p               = stats{2,5};
% % 
% % % Compute significance:
% % fprintf('Significance hit rate reduction difference (late silencing):\n  %d fast sessions vs %d slow sessions: %1.5f\n\n',sum(~isnan(FastSes)),sum(~isnan(SlowSes)),p)
% % fprintf('F(%d,%2.0f)=%1.2f, p=%1.2f; \n',stats{2,3},stats{2,4},stats{2,2},stats{2,5})
% % if p<0.05
% %     sigstar([1 2],p)
% % end
