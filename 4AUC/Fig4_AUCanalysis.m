%% This script analyzes the results from the AUROC analysis
% This approach tries to decode stimulus or decision variables based on the
% firing rate of individual neurons at certain timepoints. The analysis was performed  
% at a computer cluster. This scripts loads the result,
% displays examples, and walks through the main figures
% Oude Lohuis et al. 2022 Nat Comms

%% Reset all
startover

%% Load the data:
load('E:\Matlab\MOL_Analysis\AOudeLohuisetal_2022_NatComms\4AUC\Data4_1.mat')

%% Parameter settings:
params.Experiments          = {'ChangeDetectionConflictDecor' 'VisOnlyTwolevels' 'ChangeDetectionConflict' }; %Which versions of the task to load data from
params.ExperimentLabels     = {'NE' 'UST' 'MST'}; %Labels for each task version
params.nExperiments         = length(params.Experiments);

params                      = MOL_getColors_CHDET(params);

params.AUCbootstrapAlpha    = 0.05;

params.idx_Time             = params.xtime_AUC>0 & params.xtime_AUC<1e6;
params.minAUCBins           = 5;            %number of bins that have to be significiant


%% 





%% Figure with overlap in significant coding in V1:
params.labels_venn          = {'Ori' 'Vis' 'Hit' 'Ori-Vis'  'Ori-Hit'  'Vis-Hit' 'All'};
params.AUC_varselec         = [2 28 10];

%Make figure:
figure; set(gcf,'units','normalized','Position',[0.05 + 0.03 0.3 0.83 0.34],'color','w')

for iExp = 1:params.nExperiments
    subplot(1,3,iExp);
    %Get the index of neurons recorded in this experiment:
    idx_exp         = ismember(spikeData.session_ID,sessionData.session_ID(strcmp(sessionData.Experiment,params.Experiments(iExp))));
    idx_exp         = idx_exp & strcmp(spikeData.area,params.area);
    
    idx_sign        = []; %reinit var to store index
    frac_sign       = NaN(7,1); %store the fraction of neurons significantly coding for this variable
    %Loop over variables and find neurons that have at least n sign bins of coding for this var
    for iVar = 1:length(params.AUC_varselec)
        idx_sign(iVar,:) = sum(outputmat_sign(params.AUC_varselec(iVar),idx_exp,params.idx_Time),3)>=params.minAUCBins; %#ok<SAGROW>
    end
    
    %compute fraction of overlap for each combination:
    frac_sign(1)    = sum(idx_sign(1,:) & ~idx_sign(2,:) & ~idx_sign(3,:)) / sum(any(idx_sign));
    frac_sign(2)    = sum(~idx_sign(1,:) & idx_sign(2,:) & ~idx_sign(3,:)) / sum(any(idx_sign));
    frac_sign(3)    = sum(~idx_sign(1,:) & ~idx_sign(2,:) & idx_sign(3,:)) / sum(any(idx_sign));
    
    frac_sign(4)    = sum(idx_sign(1,:) & idx_sign(2,:) & ~idx_sign(3,:)) / sum(any(idx_sign));
    frac_sign(5)    = sum(idx_sign(1,:) & ~idx_sign(2,:) & idx_sign(3,:)) / sum(any(idx_sign));
    frac_sign(6)    = sum(~idx_sign(1,:) & idx_sign(2,:) & idx_sign(3,:)) / sum(any(idx_sign));
    frac_sign(7)    = sum(idx_sign(1,:) & idx_sign(2,:) & idx_sign(3,:)) / sum(any(idx_sign)); %triple combo
    
    frac_noresp     = sum(~any(idx_sign,1)) / sum(idx_exp);
  
    if sum(frac_sign)~=1
        error('Does not match to 100')
    end
    frac_sign(frac_sign==0)=0.01;
    
    %     [H,S] = venn(frac_sign,'FaceColor',{'r','y','b'},'FaceAlpha', 0.5,'EdgeColor','black');
    [H,S] = venn(frac_sign,'FaceColor',params.colors_AUC(params.AUC_varselec),'FaceAlpha',0.3,'EdgeColor',params.colors_AUC(params.AUC_varselec));
    frac_sign(frac_sign==0.01)=0;
    
    %Now label each zone:
    for i = 1:7
        textstring = sprintf('%s %2.0f%%',params.labels_venn{i},frac_sign(i)*100);
        text(S.ZoneCentroid(i,1), S.ZoneCentroid(i,2),textstring)
    end
    textstring = sprintf('%s %2.0f%%','Non-responsive',frac_noresp*100);
    text(-0.7,0.5,textstring)
    
    fprintf('\n(%s) Percentage of neurons coding overall:\n\n',params.ExperimentLabels{iExp})
    for i = 1:3
        fprintf('%s: %2.1f%%\n',params.labels_venn{i},(sum(idx_sign(i,:)) / sum(any(idx_sign))) *100)
    end
    
    %Print results in text
    fprintf('\nExperiment %s (total %d neurons):\n\n',params.ExperimentLabels{iExp},sum(idx_exp))
    for i = 1:7
        fprintf('%s: %3.0f neurons, %2.1f%%\n',params.labels_venn{i},frac_sign(i)* sum(idx_exp),frac_sign(i)*100)
    end
    fprintf('%s: %3.0f neurons, %2.1f%%\n','Non-responsive',frac_noresp * sum(idx_exp),frac_noresp*100)
end

%% 









%% Compute AUC over time for selected variables:
params.AUC_varselec     = [2 1 6 5 10 9];

nVars                   = length(params.AUC_varselec);
baselinetimeidx         = params.xtime_AUC<-0.05e6;

params.nTimebins_AUC = length(params.xtime_AUC);

params.auroc_fracthr    = 0.1;

frac_sign_tot           = NaN(nVars,params.nExperiments,params.nTimebins_AUC);
frac_sign_tot_norm      = NaN(nVars,params.nExperiments,params.nTimebins_AUC);

for i = 1:nVars
    iVar = params.AUC_varselec(i);
    for iExp = 1:length(params.Experiments)
        idx_exp             = ismember(spikeData.session_ID,sessionData.session_ID(strcmp(sessionData.Experiment,params.Experiments(iExp))));

        nSign               = squeeze(sum(outputmat_sign(iVar,idx_exp,:),2))';
        nTotal              = sum(~isnan(squeeze(outputmat(iVar,idx_exp,:))));
        fraction_sign       = nSign./nTotal;
        
        %store:
        frac_sign_tot(i,iExp,:) = fraction_sign;
        %normalize to max:
        fraction_sign_norm = fraction_sign - mean(fraction_sign(baselinetimeidx));
        if max(fraction_sign_norm)>params.auroc_fracthr %do not normalize if there is no coding
            fraction_sign_norm = fraction_sign_norm / max(fraction_sign_norm);
        end
        fraction_sign_norm(fraction_sign_norm<0) = 0;
        frac_sign_tot_norm(i,iExp,:) = fraction_sign_norm;
    end
end

%% Figure AUC over time (raw):

figure; set(gcf,'units','normalized','Position',[0.05 0.5 0.49 0.3],'color','w');
order = [1 4 2 5 3 6];
nVars               = length(params.AUC_varselec);
for i = 1:nVars
    iVar = params.AUC_varselec(i);
    subplot(2,nVars/2,order(i)); hold all;
    title(params.labels_AUC{iVar},'FontSize',15)
    for iExp = 1:length(params.Experiments)
        plot(params.xtime_AUC*1e-6,squeeze(frac_sign_tot(i,iExp,:)),'-','color',params.colors_experiments{iExp},'LineWidth',3);
    end
    
    if i==1
        ylabel('Fraction sign neurons','FontSize',20)
    end
    set(gca,'YTick',[0 0.5])
    ylim([0 0.5]);
    xlim([-0.5 1.5]);
%     legend(params.ExperimentLabels,'Fontsize',15); legend boxoff
end

%% Figure AUC over time (normalized):
subselec            = [2 6 10 9];
nVars               = length(subselec);

figure; set(gcf,'units','normalized','Position',[0.05 0.5 0.9 0.3],'color','w');
for i = 1:nVars
    idx = params.AUC_varselec==subselec(i);
    iVar = params.AUC_varselec(params.AUC_varselec==subselec(i));
    subplot(1,nVars,i); hold all;
    title(params.labels_AUC{iVar},'FontSize',15)
    for iExp = 1:length(params.Experiments)
        plot(params.xtime_AUC*1e-6,squeeze(frac_sign_tot_norm(idx,iExp,:)),'-','color',params.colors_experiments{iExp},'LineWidth',3);
    end
    plot([-2e6 2e6],[1 1],'--','Color',[0.5 0.5 0.5],'LineWidth',2)
    if i==1
        ylabel('Fraction sign neurons','FontSize',20)
    end
    set(gca,'YTick',[0 0.5 1])
    ylim([0 1]);
    xlim([-0.5 1.5]);
    legend(params.ExperimentLabels,'Fontsize',15); legend boxoff
end

%%











%% Make SUA heatmap over laminar depth: (visual trials)
%Histogram with binning on depth:
binedges                = 0:50:1150;
nBins                   = length(binedges)-1;

params.AUC_varselec     = [2 6 9 10]; %Get relevant coding variables
params.colormap         = {'copper' 'bone' 'pink' 'pink'}; %assign colors:

% Depth parameters:
VisGranularUpper        = 400; %Upper boundary of layer IV in V1
VisGranularLower        = 550; %Lower boundary of layer IV in V1
% Layer5Depth             = 650; %Depth from dura to center layer 5:
params.labels_layers    = {'SG' 'G' 'IG'};

for iVar = 1:length(params.AUC_varselec)
    laminardepthfig = figure; set(gcf,'units','normalized','Position',[0.01+0.24*(iVar-1) 0.5 0.22 0.32],'color','w');
    hold all; set(laminardepthfig, 'DefaultLineLineWidth', 2);
    title(params.labels_AUC{params.AUC_varselec(iVar)},'FontSize',20);
    
    idx_exp             = ismember(spikeData.session_ID,sessionData.session_ID(ismember(sessionData.Experiment,params.Experiments([2 3]))));

    AUCdepthmap         = zeros(nBins,params.nTimebins); %Init variable
    
    for iBin = 1:nBins
        idx                     = idx_exp & spikeData.ChannelY>=binedges(iBin) & spikeData.ChannelY<binedges(iBin+1);
        idx_sign                = sum(outputmat_sign(params.AUC_varselec(iVar),idx,:),2);
        AUCdepthmap(iBin,:)     = idx_sign / (sum(idx)+1);
    end
    
    AUCdepthmap(isnan(AUCdepthmap)) = 0;
    
    % Filter temporally and spatially:
    spat_filter     = fspecial('gaussian',[50 50],1.3);
    AUCdepthmap     = conv2(AUCdepthmap,spat_filter,'same');
    imagesc(flipud(AUCdepthmap));
    caxis([0 0.3])
    colormap(params.colormap{iVar});
    colorbar();

    ylabel('Depth from dura (um)','FontSize', 15)
    xlabel('Time (ms)','FontSize',15)
    
    % binedges = binedges(1:end-1) + (binedges(2)-binedges(1))/2;
    binticks = [900 550 200];
    idx = find(ismember(fliplr(binedges(1:end-1)),binticks));
    set(gca,'YTick',idx,'YTickLabel',binticks,'fontsize',15,'FontWeight','bold')
    set(gca,'fontsize',15,'FontWeight','bold')
%     timeticks = [-0.15e6 0 0.25e6 0.5e6];
    timeticks = [-0.3e6 0 0.2e6 0.5e6 1e6];
    idx = find(ismember(params.xtime,timeticks));
    set(gca,'XTick',idx,'XTickLabel',timeticks*1e-3,'fontsize',15,'FontWeight','bold')
    xlim([find(params.xtime>-0.2e6,1) find(params.xtime<0.6e6,1,'last')])
    xlim([find(params.xtime>-0.4e6,1) find(params.xtime<1.2e6,1,'last')])
%     ylim([1 nBins])
end


%% Compute for each session the fraction of significant coding neurons in each laminar zone:

idx_exp                 = ismember(spikeData.session_ID,sessionData.session_ID(ismember(sessionData.Experiment,params.Experiments([2 3]))));

uSessions               = unique(sessionData.session_ID(ismember(sessionData.Experiment,params.Experiments([2 3]))));
nSessions               = length(uSessions);

params.minNneuronsLayer = 10; %only estimate fraction with minimum of 10 neurons in this layer

%For statistics:
G_mou           = reshape(repmat(sessionData.mousename(ismember(sessionData.Experiment,params.Experiments([2 3]))),3,1),nSessions*3,1);
X_Layer         = reshape(repmat(params.labels_layers,nSessions,1),nSessions*3,1);

figure; set(gcf,'units','normalized','Position',[0.15 0.15 0.6 0.19],'color','w');
for iVar = 1:length(params.AUC_varselec)
    subplot(1,length(params.AUC_varselec),iVar); hold all;
    title(params.labels_AUC{params.AUC_varselec(iVar)},'FontSize',20);
    
    datatoplot              = NaN(length(params.labels_layers),nSessions); %Reinit variable

    switch iVar
        case 1
            params.idx_Time         = params.xtime_AUC>0e6 & params.xtime_AUC<=1e6;
        case 2
            params.idx_Time         = params.xtime_AUC>0e6 & params.xtime_AUC<=0.2e6;
        otherwise
            params.idx_Time         = params.xtime_AUC>0.2e6 & params.xtime_AUC<=1e6;
    end
    
    for iSes = 1:nSessions
        idx_ses =  ismember(spikeData.session_ID,uSessions{iSes});
        
        idx = idx_exp & spikeData.ChannelY<=VisGranularUpper & idx_ses;
        nSign = sum(outputmat_sign(params.AUC_varselec(iVar),idx,params.idx_Time),2);
        if sum(idx)>=params.minNneuronsLayer  && any(nSign~=0,3)
            datatoplot(1,iSes) = mean(nSign / (sum(idx)));
        end
        
        idx = idx_exp & spikeData.ChannelY>=VisGranularUpper & spikeData.ChannelY<=VisGranularLower & idx_ses;
        nSign = sum(outputmat_sign(params.AUC_varselec(iVar),idx,params.idx_Time),2);
        if sum(idx)>=params.minNneuronsLayer && any(nSign~=0,3)
            datatoplot(2,iSes) = mean(nSign / (sum(idx)));
        end
        
        idx = idx_exp & spikeData.ChannelY>=VisGranularLower & idx_ses;
        nSign = sum(outputmat_sign(params.AUC_varselec(iVar),idx,params.idx_Time),2);
        if sum(idx)>=params.minNneuronsLayer  && any(nSign~=0,3)
            datatoplot(3,iSes) = mean(nSign / (sum(idx)));
        end
    end
    
    meantoplot      = nanmean(datatoplot,2);
    errortoplot     = nanstd(datatoplot,[],2) ./ sqrt(sum(~isnan(datatoplot),2));
    errorbar([1 2 3],meantoplot,errortoplot,'r.','LineWidth',3);
    
    ylim([0 0.3])
    set(gca,'fontsize',15,'FontWeight','bold')
    set(gca,'XTick',[1 2 3],'XTickLabel',params.labels_layers,'fontsize',15,'FontWeight','bold')
    view([90 90])
    xlabel('Layer','FontSize', 15)
    ylabel('Fraction')
    
    %Statistical testing:
    Y_frac          = reshape(datatoplot',nSessions*3,1);
    tbl             = table(Y_frac,X_Layer,G_mou,'VariableNames',{'Fraction','Layer','Mouse'}); %Create table for mixed model
    lme             = fitlme(tbl,'Fraction~Layer+(1|Mouse)'); %construct linear mixed effects model with fixed effect of layer on fraction of coding neurons and random intercept for different mice
        
    contrasts       = {[0 1 0] [0 -1 1] [0 0 1]}; %which contrasts to do (coding is relative to first category if only one dummy)
    contrastpos     = {[1 2] [2 3] [1 3]};
    contrastlabels = {'SG vs G' 'G vs IG' 'IG vs SG'};
    
    fprintf('Posthoc comparison:\n')
    for iC = 1:3
        [p,F,DF1,DF2] = coefTest(lme,contrasts{iC});
        fprintf('(%s):F(%d,%d)=%1.2f, p=%1.3f \n',contrastlabels{iC},DF1,DF2,F,p)
        if p<0.05
            sigstar(contrastpos{iC},p);
        end
    end
    
    tempfile = fullfile(sprintf('SourceData_Fig4c_AUC_%s_Layer.xlsx',params.labels_AUC{params.AUC_varselec(iVar)}));
    writetable(tbl,tempfile)
    
%     Old statistics:
%         fprintf('Order of testing: IG:G, IG:SG, G:SG:\n')
%     combs = [1 2; 1 3; 2 3;];
%     for i = 1:3
%         p = ranksum(datatoplot(combs(i,1),:),datatoplot(combs(i,2),:));
%         if p<0.05
%             sigstar(combs(i,:),p);
%             fprintf('\nWilcoxon rank sum test - %s - layers %s vs %s - %d vs %d sessions : %1.3f\n',params.labels_AUC{params.AUC_varselec(iVar)},params.labels_layers{combs(i,1)},params.labels_layers{combs(i,2)},...
%             sum(~isnan(datatoplot(combs(i,1),:))),sum(~isnan(datatoplot(combs(i,2),:))),p)
%         end
%     end

end

%% 











%% Compute latencies:
nBoots              = 1000;

idx_resp            = params.xtime_AUC>0 & params.xtime_AUC<1e6;
idx_baseline        = params.xtime_AUC<0e6;

params.zthr         = 2;
params.fracthr      = 0.1;

params.AUC_varselec = [28 2 10 9]; %Select also vis thr hit/miss coding
nVars               = length(params.AUC_varselec);

t_crossfracthr      = NaN(params.nExperiments,nVars,nBoots); %init output vars
t_peakfracthr       = NaN(params.nExperiments,nVars,nBoots);

params.resampling = 0; subn = 128;

for iExp = 1:params.nExperiments
    idx_exp              = ismember(spikeData.session_ID,sessionData.session_ID(strcmp(sessionData.Experiment,params.Experiments(iExp))));
    fprintf('%s: %d neurons\n',params.ExperimentLabels{iExp},sum(idx_exp))

    for iVar = 1:nVars
        idx_var             = params.AUC_varselec(iVar);
        neuronidx           = find(idx_exp & ~all(isnan(squeeze(outputmat(idx_var,:,:))),2));
        if ~isempty(neuronidx)
            
            for iBoot = 1:nBoots
                %resampling for bootstrap:
                idx                 = neuronidx(randi(numel(neuronidx),numel(neuronidx),1));
%                 if params.resampling %control analysis subsampling the
%                 same amount of neurons from MST as there are UST neurons:
%                     idx                 = neuronidx(randi(numel(neuronidx),subn,1));
%                 end
                
                nSign               = squeeze(sum(outputmat_sign(idx_var,idx,:),2))';
                nTotal              = sum(~isnan(squeeze(outputmat(idx_var,idx,:))));
                
                fraction_sign       = nSign./nTotal;
                
                fracThr_boot        = nanmean(fraction_sign(idx_baseline)) + params.zthr * nanstd(fraction_sign(idx_baseline));
                
                temppoint           = params.xtime_AUC(find(fraction_sign>fracThr_boot & fraction_sign>params.fracthr & idx_resp,1));
%                 temppoint           = params.xtime(find(fraction_sign>fracThr_boot & fraction_sign>params.fracthr,1));
                
                if ~isempty(temppoint)
                    t_crossfracthr(iExp,iVar,iBoot)   = temppoint;
                    fraction_sign(~idx_resp)          = 0;
                    [~,maxidx]                        = max(fraction_sign,[],2);
                    t_peakfracthr(iExp,iVar,iBoot)    = params.xtime_AUC(maxidx);
                end
            end
        end
    end
end

t_crossfracthr(t_crossfracthr<=0) = NaN;
t_crossfracthr(t_crossfracthr>1e6) = NaN;

%% Earliest significant latency comparison Uni vs Multisensory:
plotExps            = [1 2 3];
datatoplot          = t_crossfracthr;

figure; set(gcf,'units','normalized','Position',[0.1 0.4 0.35 0.45],'color','w'); hold all;
for iExp = 1:length(plotExps)
    idx_exp             = plotExps(iExp);
    meantoplot          = squeeze(nanmean(datatoplot(:,:,:),3));
    errortoplot         = squeeze(nanstd(datatoplot(:,:,:),[],3)); % / sqrt(sum(mouseidx));

    errorbar((1:length(params.AUC_varselec))-0.3+iExp*0.2,meantoplot(idx_exp,:),errortoplot(idx_exp,:),'.','color',params.colors_experiments{iExp},'LineWidth',3,'MarkerSize',50);
    view([90 -90])
    grid on;
end

for iVar = 1:nVars
    bTrapMeanDiff = squeeze(datatoplot(2,iVar,:) - datatoplot(3,iVar,:));
    fprintf('\n')
    
    pval = 1-sum(bTrapMeanDiff<0)/nBoots;
    fprintf('%s: %1.3f\n',params.labels_AUC{params.AUC_varselec(iVar)},pval)
    if prctile(bTrapMeanDiff,100-params.AUCbootstrapAlpha/2*100)<0
        text(iVar,25e3,'*','FontSize',50)
    end
end

% Figure make up:
set(gca,'XTick',1:length(params.AUC_varselec),'XTickLabels',params.labels_AUC(params.AUC_varselec),'FontSize', 25)
set(gca,'Position',[0.25 0.2 0.7 0.7])
ylim([0 0.45e6])
set(gca,'YTick',0:0.1e6:0.4e6,'YTickLabels',(0:0.1e6:0.4e6)*1e-3,'Fontsize',15);
ylabel('Earliest time of significant fraction neurons (ms)','Fontsize',15);
legend(params.ExperimentLabels(plotExps),'Location','northeast'); legend boxoff

tempfile = fullfile('SourceData_Fig4d_Bootstrapped_AUC_Onset_Cohorts.xlsx');
for iVar = 1:nVars
    tbl = array2table(squeeze(datatoplot(:,iVar,:))','VariableNames',{'NE' 'UST' 'MST'});
    writetable(tbl,tempfile,'Sheet',iVar)
end

%% 
% meanRT_all                  = NaN(4,1); %Store the median RT for each level of visual change (2) for both visually trained task experiments (2)
bootRT_all                  = NaN(4,1000); %Store the median RT error for each level of visual change (2) for both visually trained task experiments (2)

func = @(x)(nanmean(x));
ciRT_all = NaN(4,2);

iExp = 2; %For ust, firsth thr then max
sesid                       = sessionData.session_ID(ismember(sessionData.Experiment,params.Experiments(iExp)));
[~,temptrialData]           = MOL_getTempPerSes(sesid,sessionData,trialData);%Get the sessionData for each session individually:
idx                         = strcmp(temptrialData.trialType,'X') & temptrialData.visualOriChangeNorm==2 & temptrialData.hasphotostim~=1;
fprintf('%d Vthr UST trials\n',sum(idx))
% meanRT_all(1)               = nanmean(temptrialData.responseLatency(idx));
ciRT_all(1,:)               = bootci(1000,func,temptrialData.responseLatency(idx));

for iBoot = 1:nBoots %resampling for bootstrap:
    idx_boot                 = find(idx);
    idx_boot                 = idx_boot(randi(numel(idx_boot),numel(idx_boot),1));
    bootRT_all(1,iBoot)      = nanmean(temptrialData.responseLatency(idx_boot));
end

idx                         = strcmp(temptrialData.trialType,'X') & temptrialData.visualOriChangeNorm==3 & temptrialData.hasphotostim~=1;
fprintf('%d Vmax UST trials\n',sum(idx))
% meanRT_all(2)               = nanmean(temptrialData.responseLatency(idx));
ciRT_all(2,:)               = bootci(1000,func,temptrialData.responseLatency(idx));

for iBoot = 1:nBoots %resampling for bootstrap:
    idx_boot                 = find(idx);
    idx_boot                 = idx_boot(randi(numel(idx_boot),numel(idx_boot),1));
    bootRT_all(2,iBoot)      = nanmean(temptrialData.responseLatency(idx_boot));
end

iExp = 3; %For mst, firsth thr then max
sesid                       = sessionData.session_ID(ismember(sessionData.Experiment,params.Experiments(iExp)));
[~,temptrialData]           = MOL_getTempPerSes(sesid,sessionData,trialData);%Get the sessionData for each session individually:

idx                         = strcmp(temptrialData.trialType,'X') & temptrialData.visualOriChangeNorm==2 & temptrialData.hasphotostim~=1;
fprintf('%d Vthr MST trials\n',sum(idx))
% meanRT_all(3)               = nanmean(temptrialData.responseLatency(idx));
ciRT_all(3,:)               = bootci(1000,func,temptrialData.responseLatency(idx));

for iBoot = 1:nBoots %resampling for bootstrap:
    idx_boot                 = find(idx);
    idx_boot                 = idx_boot(randi(numel(idx_boot),numel(idx_boot),1));
    bootRT_all(3,iBoot)      = nanmean(temptrialData.responseLatency(idx_boot));
end

idx                         = strcmp(temptrialData.trialType,'X') & temptrialData.visualOriChangeNorm==3 & temptrialData.hasphotostim~=1;
fprintf('%d Vmax MST trials\n',sum(idx))
% meanRT_all(4)               = nanmean(temptrialData.responseLatency(idx));
ciRT_all(4,:)               = bootci(1000,func,temptrialData.responseLatency(idx));

for iBoot = 1:nBoots %resampling for bootstrap:
    idx_boot                 = find(idx);
    idx_boot                 = idx_boot(randi(numel(idx_boot),numel(idx_boot),1));
    bootRT_all(4,iBoot)      = nanmean(temptrialData.responseLatency(idx_boot));
end

%% Make figure 4f: relationship between onset of bootstrapped hit/miss coding and reaction time:
figure; hold all; set(gcf,'units','normalized','Position',[0.05 0.4 0.25 0.37],'color','w');
mean_t_cross            = squeeze(nanmean(t_crossfracthr([2 3],[4 3],:),3));
mean_t_cross            = reshape(mean_t_cross',1,4);
lowererror_t_cross      = squeeze(prctile(t_crossfracthr([2 3],[4 3],:),(100-67)/2,3));
lowererror_t_cross      = reshape(lowererror_t_cross',1,4) - mean_t_cross;
uppererror_t_cross      = squeeze(prctile(t_crossfracthr([2 3],[4 3],:),100-(100-67)/2,3));
uppererror_t_cross      = reshape(uppererror_t_cross',1,4) - mean_t_cross;

meanRT_all              = nanmean(bootRT_all,2);
lowererrorRT_all        = ciRT_all(:,1)-meanRT_all;
uppererrorRT_all        = ciRT_all(:,2)-meanRT_all;

scatter(meanRT_all,mean_t_cross,150,'k','filled')
errorbarxy(meanRT_all,mean_t_cross,-lowererrorRT_all,uppererrorRT_all,-lowererror_t_cross,uppererror_t_cross,{'k.', 'k', 'k'})
% errorbarxy(meanRT_all,mean_t_cross,errorRT_all,error_t_cross,{'k.', 'k', 'k'})

plot([0 1e6],[200e3 200e3],'b:','LineWidth',2)

[r,p] = corr(meanRT_all,mean_t_cross');
if p<0.05
    text(150e3,320e3,sprintf('r=%1.3f,p=%1.2e',r,p),'FontSize',15)
end

% %Statistics:
MDL = fitlm(meanRT_all,mean_t_cross,'linear');
stats = anova(MDL);
fprintf('Linear regression RT on AUC onset (4 conditions) (ANOVA):\n')
fprintf('(F(%d,%2.0f) = %1.2f, p=%1.4f, ANOVA)\n',stats{1,2},stats{2,2},stats{1,4},stats{1,5})

tempfile = fullfile('SourceData_Fig4f_AUC_Onset_RT.xlsx');
writetable(table(meanRT_all,mean_t_cross','VariableNames',{'RT' 't_onset_AUC'}),tempfile)

%Figure make up:
handles = plot(MDL);
delete(handles(1));
set(handles(2),'color','k','LineWidth',1)
set(handles(3),'color','k','LineWidth',1)
set(handles(4),'color','k','LineWidth',1)
legend(gca,'off');
title('')

text(50e3,230e3,'Photodelay','FontSize',18)

%Figure make up:
xlim([100e3 600e3])
ylim([0e3 350e3])
xlabel('Reaction Time (ms)')
ylabel('Earliest AUC (ms)')
set(gca,'XTick',200e3:200e3:600e3,'XTickLabel',200:200:600,...
    'YTick',100e3:100e3:400e3,'YTickLabel',100:100:400)

%% Bootstrapping the parameters of the linear relationship between RT and onset of hit/miss coding:
figure; hold all;
p = NaN(2,nBoots);
precession = NaN(1,nBoots);

for iBoot = 1:nBoots
    xdata = bootRT_all(:,iBoot);
    ydata = reshape(t_crossfracthr([2 3],[4 3],iBoot)',1,4);
    if sum(~isnan(ydata))==4
        MDL = fit(xdata,ydata','poly1','Lower',[-inf,-inf],'Upper',[inf,inf]);
        [p(1:2,iBoot)] = coeffvalues(MDL);
        plot([0 1e6],[0*MDL.p1 + MDL.p2 1e6*MDL.p1 + MDL.p2],'r:','LineWidth',2)
        precession(iBoot) = nanmean(xdata(3:4))*MDL.p1 + MDL.p2;
    end
end


fprintf('Offset: %3.0f ms (%3.0f to %3.0f); ',median(p(2,:))*1e-3,prctile(p(2,:),5)*1e-3,prctile(p(2,:),95)*1e-3)
fprintf('Slope: %1.2f (%1.2f to %1.2f); ',median(p(1,:)),prctile(p(1,:),5),prctile(p(1,:),95))
fprintf('Mean and 95%%CI\n\n')

MDL = fit(meanRT_all,mean_t_cross','poly1','Lower',[1,-inf],'Upper',[1,inf]);
fprintf('Precession with fixed slope at 1: %3.1f ms \n\n',-MDL.p2*1e-3)

fprintf('Precession as average around reaction time: %3.1f ms \n\n',nanmean(precession)*1e-3)


%% 








