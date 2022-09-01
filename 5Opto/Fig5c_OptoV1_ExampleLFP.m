%% Opto V1 Example LFP
% This script loads an example raw trace in V1 during photostimulation of PV through ChR2
% Shown is a high pass filtered trace with 1 second bouts of photostimulation
% Reproduces Fig 5c
% Oude Lohuis et al. 2022 Nat Comms

%% 
startover

%% Load the data:
load('E:\Matlab\MOL_Analysis\AOudeLohuisetal_2022_NatComms\5Opto\Data5_1.mat')

%% Parameters: Save figure as given by user input
OutputDir           = 'E:\Documents\PhD\Figures\Project CHDET\Optogenetics\PV_Inhibition\';

Alpha               = 0.3;


%% parameters for plotting example piece:
exampletrials           = [1 5];
example_tstart          = trialData.photostimStart(exampletrials(1))-1e6;
example_tstop           = trialData.photostimStart(exampletrials(2))+2e6;

selectedchannel         = 27+32;

%% Plot for a few trials:
figure; 
lfpData.ts = lfpData.t_start:(1e6/lfpData.fs):lfpData.t_end;
set(gcf,'units','normalized','Position',[0.2 0.1 0.5 0.7],'color','w')
selectedpiece           = lfpData.ts>example_tstart & lfpData.ts<example_tstop;
plot(lfpData.ts(selectedpiece),lfpData.hpsignal{1}(selectedpiece),'k'); hold all;
ymax    = max(lfpData.hpsignal{1}(selectedpiece))*1.1;
ymin    = min(lfpData.hpsignal{1}(selectedpiece))*1.1;

for iT = exampletrials(1):exampletrials(2)
    optostart   = trialData.photostimStart(iT);
    optoend     = trialData.photostimEnd(iT);
    X = [optostart optoend optoend optostart];
    Y = [ymin ymin ymax ymax];
    patch(X, Y,'blue','EdgeColor','blue','FaceAlpha',Alpha,'EdgeAlpha',0); hold on;
end

ylim([ymin*1.2 ymax*1.2]);
xlim([example_tstart example_tstop]);
% print(gcf,fullfile(OutputDir,strcat('Ch',num2str(selectedchannel),'_Raw_Trace','.pdf')),'-dpdf','-bestfit');
% export_fig(fullfile(OutputDir,strcat('Ch7_Raw_Trace','.pdf')),'-eps');

%% Plot close up of individual pulses (not in paper:)
% Alpha = 1;

figure;
set(gcf,'units','normalized','Position',[0.2 0.1 0.5 0.7],'color','w')

preposttime         = 0.5e6;
exampletrial        = 2;
selectedpiece       = lfpData.ts>trialData.photostimStart(exampletrial)-preposttime & lfpData.ts<trialData.photostimEnd(exampletrial)+preposttime;

plot(lfpData.ts(selectedpiece),lfpData.hpsignal{1}(selectedpiece),'k','LineWidth',0.25); hold all;
ymax = max(lfpData.hpsignal{1}(selectedpiece))*1.1;
ymin = min(lfpData.hpsignal{1}(selectedpiece))*1.1;

for i = 1:1000/50
    pulsestart = trialData.photostimStart(exampletrial)+50e3*(i-1)-8e3;
    X = [pulsestart pulsestart+10e3 pulsestart+10e3 pulsestart];
    Y = [ymin ymin ymax ymax];
    patch(X, Y,'blue','EdgeColor','blue','FaceAlpha',Alpha,'EdgeAlpha',0); hold on;
end

ylim([ymin*1.1 ymax*1.1]);
set(gca,'linewidth',3)
% print(gcf,fullfile(OutputDir,strcat('Ch',num2str(selectedchannel),'_Pulse_Closeup','.pdf')),'-dpdf','-bestfit');

%% plot MUA over time: (not in paper:)
nthresholds = 4;

idx                 = find(abs(lfpData.hpsignal{1})>nthresholds*std(lfpData.hpsignal{1}));
idx                 = idx(diff(idx)~=1);

params              = params_histresponse;
params.conv_sigma   = 0.1e6;        %sd of gaussian window for smoothing

edges               = lfpData.t_start:params.binsize:lfpData.t_end;
hist_mat            = histc(lfpData.ts(idx),edges) * 1e6/params.binsize;

N                   = params.conv_twin/params.binsize;
alpha               = ((N-1)/(params.conv_sigma/params.binsize))/2; %? = (N – 1)/(2?)
win                 = gausswin(N,alpha); %convolution with gaussian
win                 = win/sum(win); %normalized

%Smooth either the total or the individual trials:
hist_mat            = padarray(hist_mat,[0 round(length(win)/2)],'symmetric','both'); %pad the array on both sides for convolution
hist_mat            = conv(hist_mat,win,'valid'); %Take only the valid overlapping center of convolution
hist_mat            = hist_mat(1:length(edges)); %slight correction to get same size (edge vs convolution)

figure;
set(gcf,'units','normalized','Position',[0.2 0.1 0.5 0.7],'color','w')
selectedpiece           = edges>example_tstart & edges<example_tstop;

plot(edges(selectedpiece),hist_mat(selectedpiece),'k'); hold all;
ymax    = max(hist_mat(selectedpiece))*1.1;
ymin    = min(hist_mat(selectedpiece))*1.1;

for iT = exampletrials(1):exampletrials(2)
    optostart   = trialData.photostimStart(iT);
    optoend     = trialData.photostimEnd(iT);
    X = [optostart optoend optoend optostart];
    Y = [ymin ymin ymax ymax];
    patch(X, Y,'blue','EdgeColor','blue','FaceAlpha',Alpha,'EdgeAlpha',0); hold on;
end

ylim([ymin*1.2 ymax*1.2]);
xlim([example_tstart example_tstop]);
set(gca,'linewidth',3)
% print(gcf,fullfile(OutputDir,strcat('Ch',num2str(selectedchannel),'_MUA_Trace','.pdf')),'-dpdf','-bestfit');

