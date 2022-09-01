% This code generates Figure 6C and/or Supplementary Figure S8g from oude Lohuis et al., 2022.
% Trials are grouped according to Side, Saliency, Modality and Outcome and
% average neuronal responses (z-scored) are computed for each condition.
% Only responsive neurons (z-score>2 for at least 1 condition) for which
% all conditions happened are kept.
% For any questions contact Umberto Olcese at u.olcese@uva.nl
clear all

%% Insert directory where OudeLohuisetal_2022_NatComms is saved:
maindirectory = 'C:\Users\jeanl\surfdrive\Shared\Manuscript - 2nd Bump\Submission - NatComm\rev1';

%% Select data
%This code uses neuronal data in Data6_2

folder = strcat(maindirectory, '\OudeLohuisetal_2022_NatComms\6VisuotactileTask\Data6_2\MST');
mstfiles = dir(strcat(folder,'\*ndata.mat'));
folder = strcat(maindirectory, '\OudeLohuisetal_2022_NatComms\6VisuotactileTask\Data6_2\UST');
ustfiles = dir(strcat(folder,'\*ndata.mat'));

%% Parameters
SelArea = 'V1';
bslwin = [-1000 -200];
rspwin = [0 1000];
thr = 2; %threshold above which z-score is considered signi

%% Initialize
clearvars frpt FRz findem abovethr

%% Obtain z-scored activity
%MST-Max-Contra-Hit
SelConds = [3];
SelTrialOutcome = [1 0 0];
[frpt{1,1}, tempsms, FRz{1,1}, ~, findem{1,1}] = JPpsth_fig6c(mstfiles, SelArea, SelConds, SelTrialOutcome);
%MST-Max-Contra-Error
SelTrialOutcome = [0 1 0];
[frpt{1,2}, tempsms, FRz{1,2}, ~, findem{1,2}] = JPpsth_fig6c(mstfiles, SelArea, SelConds, SelTrialOutcome);
%MST-Max-Contra-Miss
SelTrialOutcome = [0 0 1];
[frpt{1,3}, tempsms, FRz{1,3}, ~, findem{1,3}] = JPpsth_fig6c(mstfiles, SelArea, SelConds, SelTrialOutcome);
%MST-Max-Ipsi-Hit
SelConds = [6];
SelTrialOutcome = [1 0 0];
[frpt{1,4}, tempsms, FRz{1,4}, ~, findem{1,4}] = JPpsth_fig6c(mstfiles, SelArea, SelConds, SelTrialOutcome);
%MST-Max-Ipsi-Error
SelTrialOutcome = [0 1 0];
[frpt{1,5}, tempsms, FRz{1,5}, ~, findem{1,5}] = JPpsth_fig6c(mstfiles, SelArea, SelConds, SelTrialOutcome);
%MST-Max-Ipsi-Miss
SelTrialOutcome = [0 0 1];
[frpt{1,6}, tempsms, FRz{1,6}, ~, findem{1,6}] = JPpsth_fig6c(mstfiles, SelArea, SelConds, SelTrialOutcome);

%Same for UST
SelConds = [3];
SelTrialOutcome = [1 0 0];
[frpt{2,1}, tempsms, FRz{2,1}, ~, findem{2,1}] = JPpsth_fig6c(ustfiles, SelArea, SelConds, SelTrialOutcome);
%UST-Max-Contra-Error
SelTrialOutcome = [0 1 0];
[frpt{2,2}, tempsms, FRz{2,2}, ~, findem{2,2}] = JPpsth_fig6c(ustfiles, SelArea, SelConds, SelTrialOutcome);
%UST-Max-Contra-Miss
SelTrialOutcome = [0 0 1];
[frpt{2,3}, tempsms, FRz{2,3}, ~, findem{2,3}] = JPpsth_fig6c(ustfiles, SelArea, SelConds, SelTrialOutcome);
%UST-Max-Ipsi-Hit
SelConds = [6];
SelTrialOutcome = [1 0 0];
[frpt{2,4}, tempsms, FRz{2,4}, ~, findem{2,4}] = JPpsth_fig6c(ustfiles, SelArea, SelConds, SelTrialOutcome);
%UST-Max-Ipsi-Error
SelTrialOutcome = [0 1 0];
[frpt{2,5}, tempsms, FRz{2,5}, ~, findem{2,5}] = JPpsth_fig6c(ustfiles, SelArea, SelConds, SelTrialOutcome);
%UST-Max-Ipsi-Miss
SelTrialOutcome = [0 0 1];
[frpt{2,6}, tempsms, FRz{2,6}, ~, findem{2,6}] = JPpsth_fig6c(ustfiles, SelArea, SelConds, SelTrialOutcome);

%% Threshold and select neurons
rspwinidx = [find(tempsms>=rspwin(1),1) find(tempsms>=rspwin(2),1)];
mst_abovethr = nan(size(FRz{1,1},1), size(FRz,2));
ust_abovethr = nan(size(FRz{2,1},1), size(FRz,2));

for j = 1:min(size(FRz,2),6)  
mst_abovethr(:,j) = sum( (FRz{1,j}(:,rspwinidx(1):rspwinidx(2)) > thr) , 2) > 0 ;
ust_abovethr(:,j) = sum( (FRz{2,j}(:,rspwinidx(1):rspwinidx(2)) > thr) , 2) > 0 ;
end

%respn
mst_respn = nansum(mst_abovethr,2) > 0;
ust_respn = nansum(ust_abovethr,2) > 0;

%% Number of trials for weighing
mstW = nan(size(FRz{1,1},1),size(FRz,2) );
ustW = nan(size(FRz{2,1},1),size(FRz,2) );
for j = 1:size(FRz,2)  
    for iSess = 1: length(frpt{1,j})
        mstW( findem{1,j}(:,2) == iSess , j) = size(frpt{1,j}{iSess},1) ;
    end
    for iSess = 1: length(frpt{2,j})
        ustW( findem{2,j}(:,2) == iSess , j) = size(frpt{2,j}{iSess},1) ;
    end
end

%% Select neurons that had at least n trials in all conds.
%At least 1 trial (nan = 0)

mst_entrn_contra    = ~isnan(nanmean(FRz{1,1},2)) & ~isnan(nanmean(FRz{1,3},2)) ;
mst_entrn_ipsi      = ~isnan(nanmean(FRz{1,4},2)) & ~isnan(nanmean(FRz{1,6},2)) ; 
mst_entrn           = ~isnan(nanmean(FRz{1,1},2)) & ~isnan(nanmean(FRz{1,3},2)) ...
                    & ~isnan(nanmean(FRz{1,4},2)) & ~isnan(nanmean(FRz{1,6},2)) ; 

ust_entrn_contra    = ~isnan(nanmean(FRz{2,1},2)) & ~isnan(nanmean(FRz{2,3},2)) ;
ust_entrn_ipsi      = ~isnan(nanmean(FRz{2,4},2)) & ~isnan(nanmean(FRz{2,6},2)) ; 
ust_entrn           = ~isnan(nanmean(FRz{2,1},2)) & ~isnan(nanmean(FRz{2,3},2)) ...
                    & ~isnan(nanmean(FRz{2,4},2)) & ~isnan(nanmean(FRz{2,6},2)) ; 


%% Bootstrap for errorbars
MSTFILT = mst_respn & mst_entrn ;
USTFILT = ust_respn & ust_entrn ;
bslwinidx = [find(tempsms>=bslwin(1),1) find(tempsms>=bslwin(2),1)];
xx = [-500 1000];
subsi = [find(tempsms>=xx(1),1) find(tempsms>=xx(2),1)];

A = FRz{1,1}(MSTFILT,subsi(1):subsi(2))-nanmean(nanmean(FRz{1,1}(MSTFILT,bslwinidx(1):bslwinidx(2))));
nboot = 500;
[bci,bmeans] = bootci(nboot, {@nanmean, A}) ; %, 'Weights', mstW(mst_respn,1)
figure();
shadedErrorBar(tempsms(subsi(1):subsi(2)),nanmean(bmeans),[abs(nanmean(bmeans)-bci(2,:));abs(nanmean(bmeans)-bci(1,:))],'lineProps',{'Color',[33 7 134]/255,'LineWidth',2},'patchSaturation',0.1)
hold on

% A = FRz{1,2}(mst_respn,subsi(1):subsi(2))-nanmean(nanmean(FRz{1,2}(mst_respn,bslwinidx(1):bslwinidx(2))));
% [bci,bmeans] = bootci(nboot, {@nanmean, A}, 'alpha', 0.05, 'type','per') ;
% shadedErrorBar(tempsms(subsi(1):subsi(2)),nanmean(bmeans),[abs(nanmean(bmeans)-bci(2,:));abs(nanmean(bmeans)-bci(1,:))],'lineProps',{'Color',[.6 0 0],'LineWidth',2},'patchSaturation',0.1)

A = FRz{1,3}(MSTFILT,subsi(1):subsi(2))-nanmean(nanmean(FRz{1,3}(MSTFILT,bslwinidx(1):bslwinidx(2))));
[bci,bmeans] = bootci(nboot, {@nanmean, A}) ; %
shadedErrorBar(tempsms(subsi(1):subsi(2)),nanmean(bmeans),[abs(nanmean(bmeans)-bci(2,:));abs(nanmean(bmeans)-bci(1,:))],'lineProps',{'Color',[100 170 190]/255,'LineWidth',2},'patchSaturation',0.1)
xlim([-300 1000])

% IPSI
figure()
A = FRz{1,4}(MSTFILT,subsi(1):subsi(2))-nanmean(nanmean(FRz{1,4}(MSTFILT,bslwinidx(1):bslwinidx(2))));
[bci,bmeans] = bootci(nboot, {@nanmean, A}) ; %, 'Weights', mstW(mst_respn,4)
shadedErrorBar(tempsms(subsi(1):subsi(2)),nanmean(bmeans),[abs(nanmean(bmeans)-bci(2,:));abs(nanmean(bmeans)-bci(1,:))],'lineProps',{'Color',[33 7 134]/255,'LineWidth',2},'patchSaturation',0.1)

% A = FRz{1,5}(mst_respn,subsi(1):subsi(2))-nanmean(nanmean(FRz{1,5}(mst_respn,bslwinidx(1):bslwinidx(2))));
% [bci,bmeans] = bootci(nboot, {@nanmean, A}, 'alpha', 0.05, 'type','per') ;
% shadedErrorBar(tempsms(subsi(1):subsi(2)),nanmean(bmeans),[abs(nanmean(bmeans)-bci(2,:));abs(nanmean(bmeans)-bci(1,:))],'lineProps',{'Color',[.6 0 0],'LineWidth',2},'patchSaturation',0.1)

A = FRz{1,6}(MSTFILT,subsi(1):subsi(2))-nanmean(nanmean(FRz{1,6}(MSTFILT,bslwinidx(1):bslwinidx(2))));
[bci,bmeans] = bootci(nboot, {@nanmean, A}) ; %, 'Weights', mstW(mst_respn,6)
shadedErrorBar(tempsms(subsi(1):subsi(2)),nanmean(bmeans),[abs(nanmean(bmeans)-bci(2,:));abs(nanmean(bmeans)-bci(1,:))],'lineProps',{'Color',[100 170 190]/255,'LineWidth',2},'patchSaturation',0.1)
xlim([-300 1000])


%% Add THR
% %MST-THR-Contra-Hit
% SelConds = [2];
% SelTrialOutcome = [1 0 0];
% [frpt{3,1}, tempsms, FRz{3,1}, ~, findem{3,1}] = JPpsth_fig6c(mstfiles, SelArea, SelConds, SelTrialOutcome);
% %MST-THR-Contra-Miss
% SelTrialOutcome = [0 0 1];
% [frpt{3,3}, tempsms, FRz{3,3}, ~, findem{3,3}] = JPpsth_fig6c(mstfiles, SelArea, SelConds, SelTrialOutcome);
% %MST-THR-Ipsi-Hit
% SelConds = [5];
% SelTrialOutcome = [1 0 0];
% [frpt{3,4}, tempsms, FRz{3,4}, ~, findem{3,4}] = JPpsth_fig6c(mstfiles, SelArea, SelConds, SelTrialOutcome);
% %MST-THR-Ipsi-Miss
% SelTrialOutcome = [0 0 1];
% [frpt{3,6}, tempsms, FRz{3,6}, ~, findem{3,6}] = JPpsth_fig6c(mstfiles, SelArea, SelConds, SelTrialOutcome);
% 
% %Same for UST
% SelConds = [2];
% SelTrialOutcome = [1 0 0];
% [frpt{4,1}, tempsms, FRz{4,1}, ~, findem{4,1}] = JPpsth_fig6c(ustfiles, SelArea, SelConds, SelTrialOutcome);
% %MST-THR-Contra-Miss
% SelTrialOutcome = [0 0 1];
% [frpt{4,3}, tempsms, FRz{4,3}, ~, findem{4,3}] = JPpsth_fig6c(ustfiles, SelArea, SelConds, SelTrialOutcome);
% %MST-THR-Ipsi-Hit
% SelConds = [5];
% SelTrialOutcome = [1 0 0];
% [frpt{4,4}, tempsms, FRz{4,4}, ~, findem{4,4}] = JPpsth_fig6c(ustfiles, SelArea, SelConds, SelTrialOutcome);
% %MST-THR-Ipsi-Miss
% SelTrialOutcome = [0 0 1];
% [frpt{4,6}, tempsms, FRz{4,6}, ~, findem{4,6}] = JPpsth_fig6c(ustfiles, SelArea, SelConds, SelTrialOutcome);
% 
% %% Bootstrap for errorbars
% bslwinidx = [find(tempsms>=bslwin(1),1) find(tempsms>=bslwin(2),1)];
% xx = [-500 1000];
% subsi = [find(tempsms>=xx(1),1) find(tempsms>=xx(2),1)];
% 
% A = FRz{3,1}(mst_respn,subsi(1):subsi(2))-nanmean(nanmean(FRz{3,1}(mst_respn,bslwinidx(1):bslwinidx(2))));
% nboot = 500;
% [bci,bmeans] = bootci(nboot, {@nanmean, A}) ; %, 'Weights', mstW(mst_respn,1)
% figure();
% % errorbar(tempsms(subsi(1):subsi(2)), nanmean(bmeans), abs(nanmean(bmeans)-bci(1,:)), abs(nanmean(bmeans)-bci(2,:)))
% shadedErrorBar(tempsms(subsi(1):subsi(2)),nanmean(bmeans),[abs(nanmean(bmeans)-bci(2,:));abs(nanmean(bmeans)-bci(1,:))],'lineProps',{'Color',[0 .9 0],'LineWidth',2},'patchSaturation',0.1)
% hold on
% 
% % A = FRz{1,2}(mst_respn,subsi(1):subsi(2))-nanmean(nanmean(FRz{1,2}(mst_respn,bslwinidx(1):bslwinidx(2))));
% % [bci,bmeans] = bootci(nboot, {@nanmean, A}, 'alpha', 0.05, 'type','per') ;
% % shadedErrorBar(tempsms(subsi(1):subsi(2)),nanmean(bmeans),[abs(nanmean(bmeans)-bci(2,:));abs(nanmean(bmeans)-bci(1,:))],'lineProps',{'Color',[.6 0 0],'LineWidth',2},'patchSaturation',0.1)
% 
% A = FRz{3,3}(mst_respn,subsi(1):subsi(2))-nanmean(nanmean(FRz{3,3}(mst_respn,bslwinidx(1):bslwinidx(2))));
% [bci,bmeans] = bootci(nboot, {@nanmean, A}, 'Weights', mstW(mst_respn,3)) ; %, 'Weights', mstW(mst_respn,3)
% shadedErrorBar(tempsms(subsi(1):subsi(2)),nanmean(bmeans),[abs(nanmean(bmeans)-bci(2,:));abs(nanmean(bmeans)-bci(1,:))],'lineProps',{'Color',[.6 .6 .6],'LineWidth',2},'patchSaturation',0.1)
% 
% % IPSI
% figure()
% A = FRz{3,4}(mst_respn,subsi(1):subsi(2))-nanmean(nanmean(FRz{3,4}(mst_respn,bslwinidx(1):bslwinidx(2))));
% [bci,bmeans] = bootci(nboot, {@nanmean, A}) ; %, 'Weights', mstW(mst_respn,4)
% shadedErrorBar(tempsms(subsi(1):subsi(2)),nanmean(bmeans),[abs(nanmean(bmeans)-bci(2,:));abs(nanmean(bmeans)-bci(1,:))],'lineProps',{'Color',[0 .9 0],'LineWidth',2},'patchSaturation',0.1)
% 
% % A = FRz{1,5}(mst_respn,subsi(1):subsi(2))-nanmean(nanmean(FRz{1,5}(mst_respn,bslwinidx(1):bslwinidx(2))));
% % [bci,bmeans] = bootci(nboot, {@nanmean, A}, 'alpha', 0.05, 'type','per') ;
% % shadedErrorBar(tempsms(subsi(1):subsi(2)),nanmean(bmeans),[abs(nanmean(bmeans)-bci(2,:));abs(nanmean(bmeans)-bci(1,:))],'lineProps',{'Color',[.6 0 0],'LineWidth',2},'patchSaturation',0.1)
% 
% A = FRz{3,6}(mst_respn,subsi(1):subsi(2))-nanmean(nanmean(FRz{3,6}(mst_respn,bslwinidx(1):bslwinidx(2))));
% [bci,bmeans] = bootci(nboot, {@nanmean, A}, 'Weights', mstW(mst_respn,6)) ; %, 'Weights', mstW(mst_respn,6)
% shadedErrorBar(tempsms(subsi(1):subsi(2)),nanmean(bmeans),[abs(nanmean(bmeans)-bci(2,:));abs(nanmean(bmeans)-bci(1,:))],'lineProps',{'Color',[.6 .6 .6],'LineWidth',2},'patchSaturation',0.1)

%% Add CR and FA to Fig4c_v8
%Correct Rejection
SelConds = [1];
SelTrialOutcome = [1 1 1];
SelMouseChoice = [0 0 1]; %LRN
[frpt{1,7}, tempsms, FRz{1,7}, ~, findem{1,7}] = JPpsth_fig6c(mstfiles, SelArea, SelConds, SelTrialOutcome,SelMouseChoice);
%False Alarm Contra
SelMouseChoice = [0 1 0]; %LRN
[frpt{1,8}, tempsms, FRz{1,8}, ~, findem{1,8}] = JPpsth_fig6c(mstfiles, SelArea, SelConds, SelTrialOutcome,SelMouseChoice);
%False Alarm Ipsi
SelMouseChoice = [1 0 0]; %LRN
[frpt{1,9}, tempsms, FRz{1,9}, ~, findem{1,9}] = JPpsth_fig6c(mstfiles, SelArea, SelConds, SelTrialOutcome,SelMouseChoice);

%% PLOT EM NICELY - Fig Contra
A = FRz{1,7}(MSTFILT,subsi(1):subsi(2))-nanmean(nanmean(FRz{1,7}(MSTFILT,bslwinidx(1):bslwinidx(2))));
[bci,bmeans] = bootci(nboot, {@nanmean, A}) ; %
shadedErrorBar(tempsms(subsi(1):subsi(2)),nanmean(bmeans),[abs(nanmean(bmeans)-bci(2,:));abs(nanmean(bmeans)-bci(1,:))],'lineProps',{'Color',[148 134 15]/255,'LineWidth',2},'patchSaturation',0.1)

A = FRz{1,8}(MSTFILT,subsi(1):subsi(2))-nanmean(nanmean(FRz{1,8}(MSTFILT,bslwinidx(1):bslwinidx(2))));
[bci,bmeans] = bootci(nboot, {@nanmean, A}) ; %
shadedErrorBar(tempsms(subsi(1):subsi(2)),nanmean(bmeans),[abs(nanmean(bmeans)-bci(2,:));abs(nanmean(bmeans)-bci(1,:))],'lineProps',{'Color',[148 74 15]/255,'LineWidth',2},'patchSaturation',0.1)

%% PLOT EM NICELY - Fig Ipsi
A = FRz{1,7}(MSTFILT,subsi(1):subsi(2))-nanmean(nanmean(FRz{1,7}(MSTFILT,bslwinidx(1):bslwinidx(2))));
[bci,bmeans] = bootci(nboot, {@nanmean, A}) ; %
shadedErrorBar(tempsms(subsi(1):subsi(2)),nanmean(bmeans),[abs(nanmean(bmeans)-bci(2,:));abs(nanmean(bmeans)-bci(1,:))],'lineProps',{'Color',[148 134 15]/255,'LineWidth',2},'patchSaturation',0.1)

A = FRz{1,9}(MSTFILT,subsi(1):subsi(2))-nanmean(nanmean(FRz{1,9}(MSTFILT,bslwinidx(1):bslwinidx(2))));
[bci,bmeans] = bootci(nboot, {@nanmean, A}) ; %
shadedErrorBar(tempsms(subsi(1):subsi(2)),nanmean(bmeans),[abs(nanmean(bmeans)-bci(2,:));abs(nanmean(bmeans)-bci(1,:))],'lineProps',{'Color',[148 74 15]/255,'LineWidth',2},'patchSaturation',0.1)

%% Figure S8h (V1 responses to TACTILE stimulus. Hit/Miss)

%MST-Max-Contra-TAC
%Hit
SelConds = [8 9 10 11];
SelTrialOutcome = [1 0 0];
[frpt{1,10}, tempsms, FRz{1,10}, ~, findem{1,10}] = JPpsth_fig6c(mstfiles, SelArea, SelConds, SelTrialOutcome);
%Miss
SelConds = [8 9 10 11];
SelTrialOutcome = [0 0 1];
[frpt{1,11}, tempsms, FRz{1,11}, ~, findem{1,11}] = JPpsth_fig6c(mstfiles, SelArea, SelConds, SelTrialOutcome);

%UST-Max-Contra-TAC
%Hit
SelConds = [8 9];
SelTrialOutcome = [1 1 1];
SelMouseChoice = [0 1 0]; %LRN
[frpt{2,10}, tempsms, FRz{2,10}, ~, findem{2,10}] = JPpsth_fig6c(ustfiles, SelArea, SelConds, SelTrialOutcome);
%Miss
SelConds = [8 9];
SelTrialOutcome = [1 1 1];
SelMouseChoice = [0 0 1]; %LRN
[frpt{2,11}, tempsms, FRz{2,11}, ~, findem{2,11}] = JPpsth_fig6c(ustfiles, SelArea, SelConds, SelTrialOutcome);


%% 
figure(); hold on
A = FRz{1,10}(MSTFILT,subsi(1):subsi(2))-nanmean(nanmean(FRz{1,10}(MSTFILT,bslwinidx(1):bslwinidx(2))));
[bci,bmeans] = bootci(nboot, {@nanmean, A}) ; %
shadedErrorBar(tempsms(subsi(1):subsi(2)),nanmean(bmeans),[abs(nanmean(bmeans)-bci(2,:));abs(nanmean(bmeans)-bci(1,:))],'lineProps',{'Color',[33 7 134]/255,'LineWidth',2},'patchSaturation',0.1)

A = FRz{1,11}(MSTFILT,subsi(1):subsi(2))-nanmean(nanmean(FRz{1,11}(MSTFILT,bslwinidx(1):bslwinidx(2))));
[bci,bmeans] = bootci(nboot, {@nanmean, A}) ; %
shadedErrorBar(tempsms(subsi(1):subsi(2)),nanmean(bmeans),[abs(nanmean(bmeans)-bci(2,:));abs(nanmean(bmeans)-bci(1,:))],'lineProps',{'Color',[100 170 190]/255,'LineWidth',2},'patchSaturation',0.1)

ylim([-0.2 1])
%%
A = FRz{2,10}(USTFILT,subsi(1):subsi(2))-nanmean(nanmean(FRz{2,10}(USTFILT,bslwinidx(1):bslwinidx(2))));
[bci,bmeans] = bootci(nboot, {@nanmean, A}) ; %
shadedErrorBar(tempsms(subsi(1):subsi(2)),nanmean(bmeans),[abs(nanmean(bmeans)-bci(2,:));abs(nanmean(bmeans)-bci(1,:))],'lineProps',{'Color',[33 7 134]/255,'LineWidth',2},'patchSaturation',0.1)

A = FRz{1,11}(USTFILT,subsi(1):subsi(2))-nanmean(nanmean(FRz{2,11}(USTFILT,bslwinidx(1):bslwinidx(2))));
[bci,bmeans] = bootci(nboot, {@nanmean, A}) ; %
shadedErrorBar(tempsms(subsi(1):subsi(2)),nanmean(bmeans),[abs(nanmean(bmeans)-bci(2,:));abs(nanmean(bmeans)-bci(1,:))],'lineProps',{'Color',[100 170 190]/255,'LineWidth',2},'patchSaturation',0.1)

ylim([-0.2 1])

%% Figure S8g
%MST-VIS - Contralateral Lick (Lick Right) 
%Hit/Error/FA = VIS right / VIS left / Catch
SelLick = [0 1]; %only right licks
%Hit
SelConds = [2 3 4];
SelTrialOutcome = [1 0 0];
[frpt{4,1}, tempsms, FRz{4,1}, ~, findem{4,1}] = JPpsth_alignedtolick(mstfiles, SelArea, SelConds, SelTrialOutcome,SelLick);
%Error
SelConds = [5 6 7];
SelTrialOutcome = [0 1 0];
[frpt{4,2}, tempsms, FRz{4,2}, ~, findem{4,2}] = JPpsth_alignedtolick(mstfiles, SelArea, SelConds, SelTrialOutcome,SelLick);
%FA
SelConds = [1];
SelTrialOutcome = [1 1 1];
SelMouseChoice = [0 1 0];
[frpt{4,3}, tempsms, FRz{4,3}, ~, findem{4,3}] = JPpsth_alignedtolick(mstfiles, SelArea, SelConds, SelTrialOutcome,SelLick,SelMouseChoice);

%% 
figure(); hold on
A = FRz{4,1}(MSTFILT,subsi(1):subsi(2))-nanmean(nanmean(FRz{4,1}(MSTFILT,bslwinidx(1):bslwinidx(2))));
[bci,bmeans] = bootci(nboot, {@nanmean, A}) ; %
shadedErrorBar(tempsms(subsi(1):subsi(2)),nanmean(bmeans),[abs(nanmean(bmeans)-bci(2,:));abs(nanmean(bmeans)-bci(1,:))],'lineProps',{'Color',[33 7 134]/255,'LineWidth',2},'patchSaturation',0.1)

A = FRz{4,2}(MSTFILT,subsi(1):subsi(2))-nanmean(nanmean(FRz{4,2}(MSTFILT,bslwinidx(1):bslwinidx(2))));
[bci,bmeans] = bootci(nboot, {@nanmean, A}) ; %
shadedErrorBar(tempsms(subsi(1):subsi(2)),nanmean(bmeans),[abs(nanmean(bmeans)-bci(2,:));abs(nanmean(bmeans)-bci(1,:))],'lineProps',{'Color',[255 0 0]/255,'LineWidth',2},'patchSaturation',0.1)

A = FRz{4,3}(MSTFILT,subsi(1):subsi(2))-nanmean(nanmean(FRz{4,3}(MSTFILT,bslwinidx(1):bslwinidx(2))));
[bci,bmeans] = bootci(nboot, {@nanmean, A}) ; %
shadedErrorBar(tempsms(subsi(1):subsi(2)),nanmean(bmeans),[abs(nanmean(bmeans)-bci(2,:));abs(nanmean(bmeans)-bci(1,:))],'lineProps',{'Color',[170 170 50]/255,'LineWidth',2},'patchSaturation',0.1)

xlim([-300 200])
ylim([-0.1 1])

%% Do the same for ipsilateral condition