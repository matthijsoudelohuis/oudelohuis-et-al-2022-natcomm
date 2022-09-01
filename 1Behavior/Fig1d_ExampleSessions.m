%% ExampleSessions
% This script plots the psychometric performance of three individual example 
% sessions, one from each of the three different cohorts.
% Oude Lohuis et al. 2022 Nat Comms

%%
params.Experiments         = {'ChangeDetectionConflictDecor' 'VisOnlyPsychophysics' 'BehaviorPsychophysics'};
params.exAnimals           = {'1009' '2034' '2003'};
params.exSessions          = {'2019-03-08_14-39-42' '2019-11-16_15-28-00' '2018-01-17_13-52-00'};

%Potential other sessions for these three animals to select as example sessions:
% 1009:
% '2019-03-06_14-12-13'
% '2019-03-07_12-46-53'
% '2019-03-08_14-39-42'
% '2019-03-12_15-41-39'
% '2019-03-13_15-32-33'
% 
% 2034:
% '2019-11-14_13-45-00'
% '2019-11-15_14-04-00'
% '2019-11-16_15-28-00'
% 
% 2003:
% '2018-01-15_15-48-00'
% '2018-01-16_14-28-00'
% '2018-01-17_13-52-00'
% '2018-01-19_14-53-00'

%% Load the dataset:
load('E:\Matlab\MOL_Analysis\AOudeLohuisetal_2022_NatComms\1Behavior\Data1_1.mat','sessionData','trialData')
fprintf('Dataset: %d sessions, %d trials\n',length(sessionData.session_ID),length(trialData.session_ID));

%% Loop over animals, load the data of that session and then construct the psychometric curve
for iAnimal = 1:length(params.exAnimals)
    
    %Get the data for this session only:
    [tempsessionData,temptrialData]                 = MOL_getTempPerSes(sessionData.session_ID(iAnimal),sessionData,trialData);
    
    %use this function to get the hit rates per condition:
    [visconditions,auconditions,FullRespMat,FullnTrialsMat]     = MOL_Psy_GetTrialMat(tempsessionData,temptrialData);
    
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
    
    if numel(visconditions)==5 && numel(auconditions)==4
        FullRespMat = [FullRespMat; FullRespMat(end,:,:)]; %#ok<*AGROW>
        auconditions = [2/256 8/256 32/256 64/256 128/256];
    end
    if numel(visconditions)==6 && numel(auconditions)==4
        FullRespMat = [FullRespMat; FullRespMat(end-1:end,:,:)];
        auconditions = [1/256 2/256 8/256 32/256 64/256 128/256];
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
    % dmax = theta(1:M); n = theta(M+1:2*M); s50 = theta(2*M+1:3*M); c = theta(3*M+1:end);
    theta_est(5) = theta_est(5) * max(auconditions); %undo normalization of conditions
    theta_est(6) = theta_est(6) * max(visconditions); %undo normalization of conditions
    theta_err(5) = theta_err(5) * max(auconditions); %undo normalization of conditions
    theta_err(6) = theta_err(6) * max(visconditions); %undo normalization of conditions
    
    % Generate contingency table from fitted parameters:
    [xvals_fit_au,xvals_fit_vis,ctable_fit_mat] = MOL_Gen2ADC_PsyCurve(theta_est,params);
        
    figure; set(gcf,'color','w','units','normalized','Position', [0.1 0.1 .7 .7]);
    
    % Construct x from position of the probe trial and conditions
    xdata_au = [params.auprobepos auconditions];
    % Construct x from position of the probe trial and conditions
    xdata_vis = [params.visprobepos visconditions];
    
    set(0,'Defaultlinelinewidth',5)
    Colors_responses = {[1 0 0] [0 0 1] [0 0 0]};
    hold all;
    
    %Audio:
    subplot(1,2,1); hold all;
    %Get the hit rates for audio conditions:
    %=full first dimension, vis=1 (no change), resp=1
    datatoplot = squeeze(FullRespMat(:,1,1));
    plot(xdata_au,datatoplot,'.','Color',Colors_responses{1},'LineWidth',5,'MarkerSize',50);
    datatoplot = squeeze(FullRespMat(:,1,2));
    plot(xdata_au,datatoplot,'.','Color',Colors_responses{2},'LineWidth',5,'MarkerSize',50);
    %Plot single session fit, auditory responses first:
    plot(xvals_fit_au,squeeze(ctable_fit_mat(1,1,:)),'k','LineWidth',3);
    plot(xvals_fit_au,squeeze(ctable_fit_mat(1,2,:)),'k:','LineWidth',3);
    
    %Visual:
    subplot(1,2,2); hold all;
    %Get the hit rates for visual conditions:
    %=au=1 (no change), full visual dimension, resp=2
    datatoplot = squeeze(FullRespMat(1,:,1));
    plot(xdata_vis,datatoplot,'.','Color',Colors_responses{1},'LineWidth',5,'MarkerSize',50);
    datatoplot = squeeze(FullRespMat(1,:,2));
    plot(xdata_vis,datatoplot,'.','Color',Colors_responses{2},'LineWidth',5,'MarkerSize',50);
    %Plot single session fit, visual responses first:
    plot(xvals_fit_vis,squeeze(ctable_fit_mat(2,2,:)),'k','LineWidth',3);
    plot(xvals_fit_vis,squeeze(ctable_fit_mat(2,1,:)),'k:','LineWidth',3);
    title(sprintf('n=%d trials',length(temptrialData.trialNum)))
    MOL_Psy2Sided_FigMakeup(params)
end
