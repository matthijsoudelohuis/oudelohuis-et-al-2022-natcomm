%% This script analyzes the results from the AUROC analysis for hit/miss coding on population averaged activity
% This approach tries to decode for individual sessions the onset of hit/miss 
% coding so to relate it to RT. This scripts loads the result of the AUC analysis and reproduces Fig. 4e
% Oude Lohuis et al. 2022 Nat Comms


%% Reset all
startover

%% 
load('E:\Matlab\MOL_Analysis\AOudeLohuisetal_2022_NatComms\4AUC\Data4_2.mat','spikeData','trialData','sessionData','outputmat','outputmat','outputmat_shuf','outputmat_n')

%% Parameter settings:
params                      = params_histresponse_coding; %Parameters for PSTH (All time is in microseconds)

params.Experiments          = {'ChangeDetectionConflictDecor' 'VisOnlyTwolevels' 'ChangeDetectionConflict' }; %Which versions of the task to load data from
params.ExperimentLabels     = {'NE' 'UST' 'MST'}; %Labels for each task version
params.nExperiments         = length(params.Experiments);

params.AlignOn              = 'stimChange';      %On which timestamp to align as t=0
params.area                 = 'V1'; %Filter only V1 data

params                      = MOL_getColors_CHDET(params);

params.nshuffle             = 1000;
params.alpha                = 0.05;
params.minTrialCond         = 10;

params.minNneurons          = 10;
params.markersize           = 100;

%Firing rate parameters
params.zscore               = 1;
params.smoothing            = 1;
params.binsize              = 10e3;             %Size of the bins
params.conv_win             = 'chg';       %Type of window used for smoothing {flat, gaussian)
params.conv_sigma           = 50e3;           %sd of gaussian window for smoothing

%Construct bin edges and time axis
params.edges                = [params.t_pre:params.binsize:params.t_post] - params.binsize/2;                    %#ok<NBRAK> Define edges of the bins
params.xtime                = params.t_pre:params.binsize:params.t_post-params.binsize;                        %Define time axis
params.nTimebins            = length(params.xtime); %number of time bins

%Some measures of dimensionality:
nContrasts                  = 2;
nNeurons                    = length(spikeData.ts);
uSessions                   = unique(sessionData.session_ID);
nSessions                   = length(uSessions);

%% Construct significance:
outputmat_sig       = outputmat<prctile(outputmat_shuf,params.alpha*100,4); %Convert to significance matrix

%% Define first timebin to be significant:
t_cross_thr         = NaN(nSessions,nContrasts);

idx_resp            = params.xtime>0 & params.xtime<0.7e6;
idx_bsl             = params.xtime<0e6;

for iSes = 1:nSessions
    for iContr = 1:nContrasts
        temppoint           = params.xtime(find(outputmat_sig(iSes,:,iContr) & idx_resp,1));
       
        if ~isempty(temppoint)
            t_cross_thr(iSes,iContr)   = temppoint;
        end
    end
end

%% Compute median reaction time per session:

medianRT_all                = NaN(nSessions,nContrasts);
params.minNtrialcond        = 20;

for iSes = 1:nSessions
    sesid                       = sessionData.session_ID(iSes);
    [~,temptrialData]           = MOL_getTempPerSes(sesid,sessionData,trialData);%Get the trialData for each session individually:
    
    splits                  = {};
    splits{1}               = ismember(temptrialData.trialType,{'X'}) & temptrialData.visualOriChangeNorm==2 & temptrialData.vecResponse==2;
    splits{2}               = ismember(temptrialData.trialType,{'X'}) & temptrialData.visualOriChangeNorm==3 & temptrialData.vecResponse==2;

    if sum(splits{1})>params.minNtrialcond
        medianRT_all(iSes,1)      = nanmean(temptrialData.responseLatency(splits{1}));
    end
    if sum(splits{2})>params.minNtrialcond
        medianRT_all(iSes,2)      = nanmean(temptrialData.responseLatency(splits{2}));
    end
end
medianRT_all(medianRT_all>600e3) = NaN; %Filter out sessions with RT below over half of the response window. 

%% Make figure:
X_RT            = [];
Y_AUC           = [];

figure; set(gcf,'units','normalized','Position',[0.2 0.3 0.15 0.2],'color','w'); hold all;

%For ust sessions: (first thr then max)
iExp                = 2;
idx_exp            = strcmp(sessionData.Experiment,params.Experiments(iExp)) & outputmat_n>params.minNneurons;
xdata              = medianRT_all(idx_exp,1);
ydata              = t_cross_thr(idx_exp,1);
idx                 = ~isnan(xdata) & ~isnan(ydata); xdata = xdata(idx); ydata = ydata(idx);

scatter(xdata,ydata,params.markersize ,[1 1 1],'filled','MarkerEdgeColor',params.colors_experiments{iExp},'LineWidth',1);
H = errorbarxy(nanmean(xdata),nanmean(ydata),nanstd(xdata)/sqrt(sum(idx)),nanstd(ydata)/sqrt(sum(idx)),{'ko','k','k'});
H.hMain.Marker = 'o'; H.hMain.MarkerSize = 10; H.hMain.MarkerFaceColor = 'w'; H.hMain.MarkerEdgeColor = params.colors_experiments{iExp};
X_RT = [X_RT; xdata]; Y_AUC = [Y_AUC; ydata]; %add data to total vector for statistics

xdata              = medianRT_all(idx_exp,2);
ydata              = t_cross_thr(idx_exp,2);
idx                 = ~isnan(xdata) & ~isnan(ydata); xdata = xdata(idx); ydata = ydata(idx);
scatter(xdata,ydata,params.markersize ,params.colors_experiments{iExp},'filled','MarkerEdgeColor',params.colors_experiments{iExp},'LineWidth',1);
H = errorbarxy(nanmean(xdata),nanmean(ydata),nanstd(xdata)/sqrt(sum(idx)),nanstd(ydata)/sqrt(sum(idx)),{'ko','k','k'});
H.hMain.Marker = '.'; H.hMain.MarkerSize = 40; H.hMain.MarkerEdgeColor = params.colors_experiments{iExp};
X_RT = [X_RT; xdata]; Y_AUC = [Y_AUC; ydata];

%For mst sessions: (first thr then max)
iExp                = 3;
idx_exp            = strcmp(sessionData.Experiment,params.Experiments(iExp)) & outputmat_n>params.minNneurons;

xdata              = medianRT_all(idx_exp,1);
ydata              = t_cross_thr(idx_exp,1);
idx                 = ~isnan(xdata) & ~isnan(ydata); xdata = xdata(idx); ydata = ydata(idx);
scatter(xdata,ydata,params.markersize,[1 1 1],'filled','MarkerEdgeColor',params.colors_experiments{iExp},'LineWidth',1);
H = errorbarxy(nanmean(xdata),nanmean(ydata),nanstd(xdata)/sqrt(sum(idx)),nanstd(ydata)/sqrt(sum(idx)),{'ko','k','k'});
H.hMain.Marker = 'o'; H.hMain.MarkerSize = 10; H.hMain.MarkerFaceColor = 'w'; H.hMain.MarkerEdgeColor = params.colors_experiments{iExp};
X_RT = [X_RT; xdata]; Y_AUC = [Y_AUC; ydata];

xdata              = medianRT_all(idx_exp,2);
ydata              = t_cross_thr(idx_exp,2);
idx                 = ~isnan(xdata) & ~isnan(ydata); xdata = xdata(idx); ydata = ydata(idx);
scatter(xdata,ydata,params.markersize ,params.colors_experiments{iExp},'filled','MarkerEdgeColor',params.colors_experiments{iExp},'LineWidth',1);
H = errorbarxy(nanmean(xdata),nanmean(ydata),nanstd(xdata)/sqrt(sum(idx)),nanstd(ydata)/sqrt(sum(idx)),{'ko','k','k'});
H.hMain.Marker = '.'; H.hMain.MarkerSize = 40; H.hMain.MarkerEdgeColor = params.colors_experiments{iExp};
X_RT = [X_RT; xdata]; Y_AUC = [Y_AUC; ydata];

% %Statistics:
idx = ~(isnan(X_RT) | isnan(Y_AUC));
MDL = fitlm(X_RT(idx),Y_AUC(idx),'linear');
stats = anova(MDL);
fprintf('Linear regression RT on AUC onset (%d sessions): (ANOVA)\n',sum(idx))
fprintf('(F(%d,%2.0f) = %1.2f, p=%1.2f, ANOVA)\n',stats{1,2},stats{2,2},stats{1,4},stats{1,5})

tempfile = fullfile('SourceData_Fig4e_AUC_Onset_RT.xlsx');
writetable(table(X_RT(idx),Y_AUC(idx),'VariableNames',{'RT' 't_onset_AUC'}),tempfile)

%Figure make up:
handles = plot(MDL);
delete(handles(1));
set(handles(2),'color','k','LineWidth',1)
set(handles(3),'color','k','LineWidth',1)
set(handles(4),'color','k','LineWidth',1)
legend(gca,'off');
title('')

plot([0 1e6],[200e3 200e3],'b:','LineWidth',2)
text(110e3,230e3,'Photodelay','FontSize',18)

xlim([100e3 600e3])
ylim([-50e3 550e3])
xlabel('Reaction Time (ms)')
ylabel('Earliest AUC (ms)')
set(gca,'XTick',0:200e3:600e3,'XTickLabel',0:200:600,...
    'YTick',100e3:200e3:500e3,'YTickLabel',100:200:500)




