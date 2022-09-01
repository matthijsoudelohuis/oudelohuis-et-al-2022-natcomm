%% Example CSD
% This script plots a CSD profile of an example session
% Oude Lohuis et al. 2022 Nat Comms

%% Parameter settings:
params                      = params_histresponse(); % All time is in microseconds
params.Experiments          = {'ChangeDetectionConflictDecor' 'VisOnlyTwolevels' 'ChangeDetectionConflict' }; %Which versions of the task to load data from

params.exSession            = '2019-03-07_12-46-53';

params.savedir              = 'E:\Documents\PhD\Figures\Project CHDET\Results - 2nd Bump\2 Neural\2nd bump laminar\';

params.colormap_CSD         = 'parula';
params.sinkhot              = 1;

params.depthcorrection      = 100;

params.t_pre                = -0.1e6; %All time in microseconds
params.t_post               = .4e6;  %All time in microseconds

params.xticks               = params.t_pre:1e5:params.t_post;
params.xticklabels          = params.xticks*1e-3;

params.fs                   = 32e3; %sampling rate

params.xtime                = params.t_pre:1e6/params.fs:params.t_post(1) - 1e6/params.fs;
params.nTimebins            = numel(params.xtime);

%% figure settings:
set(0,'defaultLineLineWidth',2)
set(0,'defaultaxesfontsize',14)
set(0,'DefaultAxesFontName','Arial')
 
%% Get data for example session:
[Data] = MOL_GetData('E:','CHDET',params.Experiments,[],params.exSession,{'sessionData' 'trialData' 'lfpData_lazy'});
sessionData     = Data.sessionData;
trialData       = Data.trialData;
lfpData         = Data.lfpData;

%% Filter channels on area:
idx             = ismember(lfpData.area,{'V1'});
lfpFields     = fieldnames(lfpData);
for iF = 1:length(lfpFields)
    lfpData.(lfpFields{iF}) = lfpData.(lfpFields{iF})(idx,:);
end

%%  MUA power:
figure; hold all; set(gcf,'units','normalized','Position',[0.12 0.45 0.41 0.34],'color','w');



%Get power estimate for different frequency bands:
[tempsessionData,temptrialData,templfpData] = MOL_getTempPerSes(sessionData.session_ID(1),sessionData,trialData,lfpData);

templfpData.sortedChannelDepth = -templfpData.sortedChannelDepth;
templfpData.sortedChannelDepth = templfpData.sortedChannelDepth + params.depthcorrection;
%Panel 1: MUA power
subplot(1,2,1);

plot(templfpData.HF_PWR(templfpData.sortedisgood==1),templfpData.sortedChannelDepth(templfpData.sortedisgood==1),'k.-','LineWidth',2,'MarkerSize',25)
xlabel('Normalized MUA power')
ylim([-1000 0])
set(gca,'YDir','normal')

%Panel 2: CSD power
subplot(1,2,2)
meancsd             = templfpData.CheckerCSD;
meanerp            = templfpData.CheckerERP;

params.cscale               = [-max(max(meancsd))*0.95 max(max(meancsd))*0.95];

imagesc(params.xtime,templfpData.sortedChannelDepth,meancsd,params.cscale); hold on;
switch params.colormap_CSD
    case 'redblue'
        h       = coolwarm + repmat(0.1-abs(linspace(-0.1,.1,64))',1,3);
        h(h>1)  = 1;
        if params.sinkhot
            h = flipud(h);
        end
        colormap(h);
    case 'parula'
        if params.sinkhot
            colormap(flipud(parula));
        else                     colormap(parula);
        end
end

offsetmat = repmat(templfpData.sortedChannelDepth,1,length(params.xtime));
% plot(params.xtime,meanerp*20e4 + offsetmat,'k','LineWidth',0.5); hold on;
idx_subsample = 1:10:params.nTimebins;
plot(params.xtime(idx_subsample),meanerp(:,idx_subsample)*10e4 + offsetmat(:,idx_subsample),'k','LineWidth',0.5); hold on;

ylabel('Channel','FontSize', 15)
plot([0 0],ylim,'k','LineWidth',1);
xlim([params.xtime(1) params.xtime(end)]);
set(gca,'XTick',params.xticks,'XTickLabels',params.xticklabels)
xlabel('Time from stimulus (ms)','FontSize', 15)

ylim([-1000 0])
set(gca,'YDir','normal')

% export_fig(fullfile(params.savedir,sprintf('ExSession_%s',tempsessionData.session_ID{1})),'-eps','-nocrop')
export_fig(fullfile(params.savedir,sprintf('CSD_MUA_ExSession_%s',tempsessionData.session_ID{1})),'-eps')
