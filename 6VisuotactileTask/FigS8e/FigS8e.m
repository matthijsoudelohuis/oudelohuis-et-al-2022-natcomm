%% This code regens Fig S8e
%% f-Psychometric curves with control, early, late silencing

% MST-Early : J016 and 019 are not OptoLR, 21 and 22 are.
folder = 'C:\Scratch\Beh\Data\V1 second bump silencing\V1 full stim silencing';
files = dir(strcat(folder,'\*.mat'));
[TotalTrials_me,HitLeft_me,HitRight_me,~,NumTr_me] = getHitRates(files);
%For non-OptoLR, remove conditions 40 and 42, since they are actually
%multisensory
HitLeft_me(1:4,40) = [nan; nan; nan; nan];
HitLeft_me(1:4,42) = [nan; nan; nan; nan];
HitRight_me(1:4,40) = [nan; nan; nan; nan];
HitRight_me(1:4,42) = [nan; nan; nan; nan];
%MST-Late : J016 and 019 are not OptoLR, 21 and 22 are.
folder = 'C:\Scratch\Beh\Data\V1 second bump silencing\V1 delayed silencing';
files = dir(strcat(folder,'\*.mat'));
[TotalTrials_ml,HitLeft_ml,HitRight_ml,RT,NumTr_ml] = getHitRates(files);
%For non-OptoLR, remove conditions 40 and 42, since they are actually
%multisensory
HitLeft_ml(1:3,40) = [nan; nan; nan];
HitLeft_ml(1:3,42) = [nan; nan; nan];
HitRight_ml(1:3,40) = [nan; nan; nan];
HitRight_ml(1:3,42) = [nan; nan; nan];
%UST-Both : all are LRD
folder = 'C:\Scratch\Beh\Data\Second Bump Unisensory - selection Jean 4';
files = dir(strcat(folder,'\*.mat'));
[TotalTrials_u,HitLeft_u,HitRight_u,~,NumTr_u] = getHitRates(files);

%% Build psychometrics (+sem)
figure(); hold on
%MST 
    %MST - Control - Hit Right - Visual  
    HR_Ctrl = [HitRight_me(:,[6 5 1 2 3]) ; HitRight_ml(:,[6 5 1 2 3])]
    errorbar([1 2 3 4 5],nanmean(HR_Ctrl), nanstd(HR_Ctrl)./sqrt(sum(~isnan(HR_Ctrl))),'o-', 'Color',[.3 0 0])
    ylim([0 100])
    %MST - Early - Hit Right - Visual
    HR_Early = [HitRight_me(:,[40 16 17 18])];
    errorbar([2 3 4 5], nanmean(HR_Early), nanstd(HR_Early)./sqrt(sum(~isnan(HR_Early))),'o-', 'Color', [1 0 0])
    %MST - Late - Hit Right - Visual
    HR_Late = [HitRight_ml(:,[40 16 17 18])];
    errorbar([2 3 4 5], nanmean(HR_Late), nanstd(HR_Late)./sqrt(sum(~isnan(HR_Late))),'o-', 'Color', [.6 0 0])

    %MST - Control - Hit LEFT - Visual  
    HR_Ctrl = [HitLeft_me(:,[6 5 1 2 3]) ; HitLeft_ml(:,[6 5 1 2 3])]
    errorbar([1 2 3 4 5],nanmean(HR_Ctrl), nanstd(HR_Ctrl)./sqrt(sum(~isnan(HR_Ctrl))),'o-', 'Color',[0 .3 0])
    ylim([0 100])
    %MST - Early - Hit LEFT - Visual
    HR_Early = [HitLeft_me(:,[40 16 17 18])];
    errorbar([2 3 4 5], nanmean(HR_Early), nanstd(HR_Early)./sqrt(sum(~isnan(HR_Early))),'o-', 'Color', [0 1 0])
    %MST - Late - Hit LEFT - Visual
    HR_Late = [HitLeft_ml(:,[40 16 17 18])];
    errorbar([2 3 4 5], nanmean(HR_Late), nanstd(HR_Late)./sqrt(sum(~isnan(HR_Late))),'o-', 'Color', [0 .6 0])

%UST
figure();
hold on
    %MST - Control - Hit Right - Visual  
    HR_Ctrl = [HitRight_u(:,[6 5 1 2 3])];
    errorbar([1 2 3 4 5],nanmean(HR_Ctrl), nanstd(HR_Ctrl)./sqrt(sum(~isnan(HR_Ctrl))),'o-', 'Color',[.3 0 0])
    ylim([0 100])
    %MST - Early - Hit Right - Visual
    HR_Early = [HitRight_u(:,[21 20 16 17 18])];
    errorbar([1 2 3 4 5], nanmean(HR_Early), nanstd(HR_Early)./sqrt(sum(~isnan(HR_Early))),'o-', 'Color', [1 0 0])
    %MST - Late - Hit Right - Visual
    HR_Late = [HitRight_u(:,[31 30 26 27 28])];
    errorbar([1 2 3 4 5], nanmean(HR_Late), nanstd(HR_Late)./sqrt(sum(~isnan(HR_Late))),'o-', 'Color', [.6 0 0])

    %MST - Control - Hit Left - Visual  
    HR_Ctrl = [HitLeft_u(:,[6 5 1 2 3])];
    errorbar([1 2 3 4 5],nanmean(HR_Ctrl), nanstd(HR_Ctrl)./sqrt(sum(~isnan(HR_Ctrl))),'o-', 'Color',[0 .3 0])
    ylim([0 100])
    %MST - Early - Hit Left - Visual
    HR_Early = [HitLeft_u(:,[21 20 16 17 18])];
    errorbar([1 2 3 4 5], nanmean(HR_Early), nanstd(HR_Early)./sqrt(sum(~isnan(HR_Early))),'o-', 'Color', [0 1 0])
    %MST - Late - Hit Left - Visual
    HR_Late = [HitLeft_u(:,[31 30 26 27 28])];
    errorbar([1 2 3 4 5], nanmean(HR_Late), nanstd(HR_Late)./sqrt(sum(~isnan(HR_Late))),'o-', 'Color', [0 .6 0])

        %% STATIES FOR %LICK RIGHT
%MST Contra Thr    
    p_thr_mst_e = signrank( HitRight_me(:,[2]), HitRight_me(:,[17]) );
    p_thr_mst_l = signrank( HitRight_ml(:,[2]), HitRight_ml(:,[17]) );
%UST Contra Thr    
    p_thr_ust_e = signrank( HitRight_u(:,[2]), HitRight_u(:,[17]) ) ;
    p_thr_ust_l = signrank( HitRight_u(:,[2]), HitRight_u(:,[27]) ) ;
    
%MST FA 
    p_far_mst_e = signrank( HitRight_me(:,[1]), HitRight_me(:,[16]) );
    p_far_mst_l = signrank( HitRight_ml(:,[1]), HitRight_ml(:,[16]) );
%UST FA 
    p_far_ust_e = signrank( HitRight_u(:,[1]), HitRight_u(:,[16]) ) ;
    p_far_ust_l = signrank( HitRight_u(:,[1]), HitRight_u(:,[26]) ) ;    
    
pp = signrank( HitRight_u(:,[17]), HitRight_u(:,[27]) )
    
% %MST Contra Max    
    p_max_mst_e = signrank( HitRight_me(:,[3]), HitRight_me(:,[18]) );
    p_max_mst_l = signrank( HitRight_ml(:,[3]), HitRight_ml(:,[18]) );
%UST Contra Max    
    p_max_ust_e = signrank( HitRight_u(:,[3]), HitRight_u(:,[18]) ) ;
    p_max_ust_l = signrank( HitRight_u(:,[3]), HitRight_u(:,[28]) ) ;

%ps4cor = [p_thr_mst_e p_thr_mst_l p_thr_ust_e p_thr_ust_l];
% % FDR
% [Q] = mafdr(ps4cor,'BHFDR','true')
% Q<0.05

% Bonferroni-Holm
corp_u1 = bonf_holm([p_thr_ust_e p_thr_ust_l  p_thr_mst_e p_thr_mst_l]);

%% h-Dprime silencing results with laser on headbar (control).
 folder = 'C:\Scratch\Beh\Data\Control laser on s1 or headbar\UST';
 files = dir(strcat(folder,'\**\*.mat'));
 [TotalTrials,HitLeft,HitRight,MedianRT,NumTr] = getHitRates(files);
% HitRight(HitRight==0) = 100*((HitRight(HitRight==0)./100).*TotalTrials(HitRight==0) +1)./(TotalTrials(HitRight==0)+2); 
% HitLeft(HitLeft==0) = 100*((HitLeft(HitLeft==0)./100).*TotalTrials(HitLeft==0) +1)./(TotalTrials(HitLeft==0)+2) ; 
% HitRight(HitRight==100) = 100*((HitRight(HitRight==100)./100).*TotalTrials(HitRight==100) +1)./(TotalTrials(HitRight==100)+2); 
% HitLeft(HitLeft==100) = 100*((HitLeft(HitLeft==100)./100).*TotalTrials(HitLeft==100) +1)./(TotalTrials(HitLeft==100)+2) ; 
figure();
hold on
    %UST - Control - Hit Right - Visual  
    HR_Ctrl = [HitRight(:,[7 6 5 1 2 3 4])];
    errorbar([1 2 3 4 5 6 7],nanmean(HR_Ctrl,1), nanstd(HR_Ctrl,[],1)./sqrt(sum(~isnan(HR_Ctrl))),'o-', 'Color',[.3 0 0])
    ylim([0 100])
    %UST - Early - Hit Right - Visual
    HR_Early = nan(8,6);
    HR_Early(7:8,:) = HitRight([4,7],[21 20 16 17 18 19]);
    HR_Early(1:6,2:end) = HitRight([1 2 3 5 6 8],[40 16 17 18 37]);
    errorbar([2 3 4 5 6 7], nanmean(HR_Early,1), nanstd(HR_Early,[],1)./sqrt(sum(~isnan(HR_Early))),'o-', 'Color', [1 0 0])
    
    %UST - Control - Hit Left - Visual  
    HL_Ctrl = [HitLeft(:,[7 6 5 1 2 3 4])];
    errorbar([1 2 3 4 5 6 7],nanmean(HL_Ctrl,1), nanstd(HL_Ctrl,[],1)./sqrt(sum(~isnan(HL_Ctrl))),'o-', 'Color',[0 .3 0])
    ylim([0 100])
    %UST - Early - Hit Left - Visual
    HL_Early = nan(8,6);
    HL_Early(7:8,:) = HitLeft([4,7],[21 20 16 17 18 19]);
    HL_Early(1:6,2:end) = HitLeft([1 2 3 5 6 8],[40 16 17 18 37]);
    errorbar([2 3 4 5 6 7], nanmean(HL_Early,1), nanstd(HL_Early,[],1)./sqrt(sum(~isnan(HL_Early))),'o-', 'Color', [0 1 0])

 %% MST
 folder = 'C:\Scratch\Beh\Data\Control laser on s1 or headbar\MST';
 files = dir(strcat(folder,'\**\*.mat'));
 [TotalTrials,HitLeft,HitRight,MedianRT,NumTr] = getHitRates(files);
% HitRight(HitRight==0) = 100*((HitRight(HitRight==0)./100).*TotalTrials(HitRight==0) +1)./(TotalTrials(HitRight==0)+2); 
% HitLeft(HitLeft==0) = 100*((HitLeft(HitLeft==0)./100).*TotalTrials(HitLeft==0) +1)./(TotalTrials(HitLeft==0)+2) ; 
% HitRight(HitRight==100) = 100*((HitRight(HitRight==100)./100).*TotalTrials(HitRight==100) +1)./(TotalTrials(HitRight==100)+2); 
% HitLeft(HitLeft==100) = 100*((HitLeft(HitLeft==100)./100).*TotalTrials(HitLeft==100) +1)./(TotalTrials(HitLeft==100)+2) ; 

 figure()
 hold on
     %MST - Control - Hit Right - Visual  
    HR_Ctrl = [HitRight(:,[6 5 1 2 3])];
    errorbar([1 2 3 4 5],nanmean(HR_Ctrl), nanstd(HR_Ctrl)./sqrt(sum(~isnan(HR_Ctrl))),'o-', 'Color',[.3 0 0])
    ylim([0 100])
    %MST - Early - Hit Right - Visual
    HR_Early = [HitRight(:,[16 17 18])];
    errorbar([ 3 4 5], nanmean(HR_Early), nanstd(HR_Early)./sqrt(sum(~isnan(HR_Early))),'o-', 'Color', [1 0 0])
    
    %MST - Control - Hit lEFT - Visual  
    HL_Ctrl = [HitLeft(:,[6 5 1 2 3])];
    errorbar([1 2 3 4 5],nanmean(HL_Ctrl), nanstd(HL_Ctrl)./sqrt(sum(~isnan(HL_Ctrl))),'o-', 'Color',[0 .3 0])
    %MST - Early - Hit lEFT - Visual
    HL_Early = [HitLeft(:,[16 17 18])];
    errorbar([ 3 4 5], nanmean(HL_Early), nanstd(HL_Early)./sqrt(sum(~isnan(HL_Early))),'o-', 'Color', [0 1 0])
    