% This code generates Supplementary Figures S8a,b from oude Lohuis et al., 2022.
% Figure S8a- Dprimes
% Use only Max conditions

%% Insert directory where OudeLohuisetal_2022_NatComms is saved:
maindirectory = 'C:\Users\jeanl\surfdrive\Shared\Manuscript - 2nd Bump\Submission - NatComm\rev2';

%% Get data
%This code obtains data from Data6_3

%UST
folder = strcat(maindirectory,'\OudeLohuisetal_2022_NatComms\6VisuotactileTask\Data6_3\UST');
files = dir(strcat(folder,'\*.mat'));
format = 'LRDual';
[DprimesR_UST, MedianRT_UST, NumTrUST, MouseID_ust] = getDprimes(files,format,1,0);
[DprimesL_UST,~,~] = getDprimes(files,format,1,1);

D_V_UST = [DprimesR_UST(:,3) ; DprimesL_UST(:,6)];
D_T_UST = [DprimesR_UST(:,9) ; DprimesL_UST(:,11)];


%MST - Early and Late
folder = strcat(maindirectory,'\OudeLohuisetal_2022_NatComms\6VisuotactileTask\Data6_3\MST\V1 delayed silencing');
files1 = dir(strcat(folder,'\*.mat'));
folder2 = strcat(maindirectory,'\OudeLohuisetal_2022_NatComms\6VisuotactileTask\Data6_3\MST\V1 full stim silencing');
files2 = dir(strcat(folder2,'\*.mat'));
files = [files1; files2];
format = 'Rec2';
[DprimesR_MST, MedianRT_MST, NumTrMST, MouseID_mst] = getDprimes(files,format,1,0);
[DprimesL_MST,~,~] = getDprimes(files,format,1,1);

D_V_MST = [DprimesR_MST(:,3) ; DprimesL_MST(:,6)];
D_T_MST = [DprimesR_MST(:,9) ; DprimesL_MST(:,11)];

%% New stats : LME

%Build table
    %UST part
    Dprime = [DprimesR_UST(:,3) ; DprimesL_UST(:,6)];
    Mouse  = [MouseID_ust; MouseID_ust];
    tbl1 = table(Dprime,Mouse);
    tbl1.Context = repmat({'UST'},size(tbl1,1),1);
    tbl1.Side = [repmat({'R'},size(DprimesR_UST(:,3),1),1); repmat({'L'},size(DprimesL_UST(:,6),1),1)] ; 
    %MST part
    Dprime = [DprimesR_MST(:,3) ; DprimesL_MST(:,6)];
    Mouse = [MouseID_mst; MouseID_mst];
    tbl2 = table(Dprime,Mouse);
    tbl2.Context = repmat({'MST'},size(tbl2,1),1);
    tbl2.Side = [repmat({'R'},size(DprimesR_MST(:,3),1),1); repmat({'L'},size(DprimesL_MST(:,6),1),1)] ; 

%Put together
tbl = [tbl1; tbl2];
tbl.Context       = categorical(tbl.Context);
tbl.Mouse         = categorical(tbl.Mouse);
tbl.Side          = categorical(tbl.Side);

%LME model
% lme = fitlme(tbl,'Dprime ~ Context + Side + (1|Mouse:Context)');
lme = fitlme(tbl,'Dprime ~ Context + Side');
stats = dataset2table(anova(lme,'DFMethod','Satterthwaite'));
p1 = stats{2,5};


%% Plot panel S8a : Dprimes
figure(); hold on

scatter(1*ones(size(D_V_UST)),D_V_UST, 'filled')
errorbar(2,nanmedian(D_V_UST),iqr(D_V_UST)/2,'o')
scatter(3*ones(size(D_T_UST)),D_T_UST, 'filled')
errorbar(4,nanmedian(D_T_UST),iqr(D_T_UST)/2,'o')

scatter(5*ones(size(D_V_MST)),D_V_MST, 'filled')
errorbar(6,nanmedian(D_V_MST),iqr(D_V_MST)/2,'o')
scatter(7*ones(size(D_T_MST)),D_T_MST, 'filled')
errorbar(8,nanmedian(D_T_MST),iqr(D_T_MST)/2,'o')

xlim([0 9])
set(gca,'TickDir','out')
ylabel('D-prime')
xlabel('UST visual - UST tactile - MST visual - MST tactile')

%% Figure S9i (removed) - D'vsRT
figure(); hold on
%UST
scatter(MedianRT_UST(:,2),DprimesR_UST(:,2),[],[0 0 .5],'filled') %Low Right UST
scatter(MedianRT_UST(:,3),DprimesR_UST(:,3),[],[0 0 1],'filled') %High Right UST
%MST
scatter(MedianRT_MST(:,2),DprimesR_MST(:,2),[],[0 0 .5],'filled') %Low Right UST
scatter(MedianRT_MST(:,3),DprimesR_MST(:,3),[],[0 0 1],'filled') %High Right UST
RT = [MedianRT_UST(:,2);MedianRT_UST(:,3);MedianRT_MST(:,2);MedianRT_MST(:,3)];
D = [DprimesR_UST(:,2);DprimesR_UST(:,3);DprimesR_MST(:,2);DprimesR_MST(:,3)];

[r,p] = corrcoef(RT,D)
xlabel('RT')
ylabel('Dprime')

%% Figure S9b - Thresholds
% Use logistic fits. Define threshold as the value between FAR and
% saturation (max %Correct)(half-way)
thresholds = struct;
modeldprimes = struct;

%% UST 
%Visual
folder = strcat(maindirectory,'\OudeLohuisetal_2022_NatComms\6VisuotactileTask\Data6_1\UST');
files = dir(strcat(folder,'\**\All\MLogisticReg_Visual*.mat'));
% Visual RIGHT
for iMuis = 1:length(files)
    load(fullfile(files(iMuis).folder,files(iMuis).name)) 
    %Need to scale to 0.01 till maximum
    xR = [0:0.01:0.99 1:0.5:100].^n;
    xL = xR;
    %Recalculate pR and pL for new x-axis
    pL=[];pR=[];RHS_L=[];RHS_R=[];
    for i = 1:length(xL)
        for j = 1:length(xR)
            RHS_L(i,j) = B(1,1) + B(2,1)*xL(i) + B(3,1)*xR(j); 
            RHS_R(i,j) = B(1,2) + B(2,2)*xL(i) + B(3,2)*xR(j); 
            pR(i,j) = exp(RHS_R(i,j)) ./ (1+ exp(RHS_L(i,j)) + exp(RHS_R(i,j)));
        end
    end
    xxR = xR.^(1/n);
    half = pR(1,1) + (pR(1,end)-pR(1,1))/2; %Half-detection value
    thresholds(iMuis).UVR = xxR(find(pR(1,:)>half,1)); %Visual Threshold Right
    modeldprimes(iMuis).UVR = norminv(pR(1,end))-norminv(pR(1,1));
end
% Visual LEFT
for iMuis = 1:length(files)
    load(fullfile(files(iMuis).folder,files(iMuis).name)) 
    %Need to scale to 0.01 till maximum
    xR = [0:0.01:0.99 1:0.5:100].^n;
    xL = xR;
    %Recalculate pR and pL for new x-axis
    pL=[];pR=[];RHS_L=[];RHS_R=[];
    for i = 1:length(xL)
        for j = 1:length(xR)
            RHS_L(i,j) = B(1,1) + B(2,1)*xL(i) + B(3,1)*xR(j); 
            RHS_R(i,j) = B(1,2) + B(2,2)*xL(i) + B(3,2)*xR(j); 
            pL(i,j) = exp(RHS_L(i,j)) ./ (1+ exp(RHS_L(i,j)) + exp(RHS_R(i,j)));
        end
    end    
    xxL = xL.^(1/n);
    half = pL(1,1) + (pL(end,1)-pL(1,1))/2; %Half-detection value
    thresholds(iMuis).UVL = xxL(find(pL(:,1)>half,1)); %Visual Threshold Left
    modeldprimes(iMuis).UVL = norminv(pL(end,1))-norminv(pL(1,1));
end
%Tactile - Makes no sense but still to see
files = dir(strcat(folder,'\**\All\MLogisticReg_Tactile*.mat'));
%Tactile RIGHT
for iMuis = 1:length(files)
    load(fullfile(files(iMuis).folder,files(iMuis).name)) 
    xR = [0:0.01:0.99 1:0.1:9.9 10:1:100].^n; xL=xR;
    %Recalculate pR and pL for new x-axis
    pL=[];pR=[];RHS_L=[];RHS_R=[];
    for i = 1:length(xL)
        for j = 1:length(xR)
            RHS_L(i,j) = B(1,1) + B(2,1)*xL(i) + B(3,1)*xR(j); 
            RHS_R(i,j) = B(1,2) + B(2,2)*xL(i) + B(3,2)*xR(j); 
            pR(i,j) = exp(RHS_R(i,j)) ./ (1+ exp(RHS_L(i,j)) + exp(RHS_R(i,j)));
        end
    end 
    xxR = xR.^(1/n);
    half = pR(1,1) + (pR(1,end)-pR(1,1))/2; %Half-detection value
    thresholds(iMuis).UTR = xxR(find(pR(1,:)>half,1)); %Visual Threshold Right
    modeldprimes(iMuis).UTR = norminv(pR(1,end))-norminv(pR(1,1));
end
%Tactile LEFT
for iMuis = 1:length(files)
    load(fullfile(files(iMuis).folder,files(iMuis).name)) 
    xR = [0:0.01:0.99 1:0.1:9.9 10:1:100].^n; xL=xR;
    %Recalculate pR and pL for new x-axis
    pL=[];pR=[];RHS_L=[];RHS_R=[];
    for i = 1:length(xL)
        for j = 1:length(xR)
            RHS_L(i,j) = B(1,1) + B(2,1)*xL(i) + B(3,1)*xR(j); 
            RHS_R(i,j) = B(1,2) + B(2,2)*xL(i) + B(3,2)*xR(j); 
            pL(i,j) = exp(RHS_L(i,j)) ./ (1+ exp(RHS_L(i,j)) + exp(RHS_R(i,j)));
        end
    end 
    xxL = xL.^(1/n);
    half = pL(1,1) + (pL(end,1)-pL(1,1))/2; %Half-detection value
    thresholds(iMuis).UTL = xxL(find(pL(:,1)>half,1)); %Visual Threshold Left
    modeldprimes(iMuis).UTL = norminv(pL(end,1))-norminv(pL(1,1));    
end
%% MST 
%Visual
folder = strcat(maindirectory,'\OudeLohuisetal_2022_NatComms\6VisuotactileTask\Data6_1\MST');
files = dir(strcat(folder,'\**\all\MLogisticReg_Visual_nomaxcut*.mat'));
%For this fit we only use the behavioral data from mice used in the
%silencing experiment:
files = files([10 13 15 16]);
% Visual RIGHT
for iMuis = 1:length(files)
    load(fullfile(files(iMuis).folder,files(iMuis).name)) 
    %Need to scale to 0.01 till maximum
    xR = [0:0.01:0.99 1:0.5:100].^n;
    xL = xR;
    %Recalculate pR and pL for new x-axis
    pL=[];pR=[];RHS_L=[];RHS_R=[];
    for i = 1:length(xL)
        for j = 1:length(xR)
            RHS_L(i,j) = B(1,1) + B(2,1)*xL(i) + B(3,1)*xR(j); 
            RHS_R(i,j) = B(1,2) + B(2,2)*xL(i) + B(3,2)*xR(j); 
            pR(i,j) = exp(RHS_R(i,j)) ./ (1+ exp(RHS_L(i,j)) + exp(RHS_R(i,j)));
        end
    end    
    xxR = xR.^(1/n); 
    half = pR(1,1) + (pR(1,end)-pR(1,1))/2; %Half-detection value
    thresholds(iMuis).MVR = xxR(find(pR(1,:)>half,1)); %Visual Threshold Right
    modeldprimes(iMuis).MVR = norminv(pR(1,end))-norminv(pR(1,1));
end
% Visual LEFT
for iMuis = 1:length(files)
    load(fullfile(files(iMuis).folder,files(iMuis).name)) 
    %Need to scale to 0.01 till maximum
    xR = [0:0.01:0.99 1:0.5:100].^n;
    xL = xR;
    %Recalculate pR and pL for new x-axis
    pL=[];pR=[];RHS_L=[];RHS_R=[];
    for i = 1:length(xL)
        for j = 1:length(xR)
            RHS_L(i,j) = B(1,1) + B(2,1)*xL(i) + B(3,1)*xR(j); 
            RHS_R(i,j) = B(1,2) + B(2,2)*xL(i) + B(3,2)*xR(j); 
            pL(i,j) = exp(RHS_L(i,j)) ./ (1+ exp(RHS_L(i,j)) + exp(RHS_R(i,j)));
        end
    end    
    xxL = xL.^(1/n); 
    half = pL(1,1) + (pL(end,1)-pL(1,1))/2; %Half-detection value
    thresholds(iMuis).MVL = xxL(find(pL(:,1)>half,1)); %Visual Threshold Left
    modeldprimes(iMuis).MVL = norminv(pL(end,1))-norminv(pL(1,1));
end
%Tactile - Makes no sense but still to see
files = dir(strcat(folder,'\**\all\MLogisticReg_Tactile*.mat'));
files = files([10 13 15 16]);
%Tactile RIGHT
for iMuis = 1:length(files)
    load(fullfile(files(iMuis).folder,files(iMuis).name)) 
    xR = [0:0.01:0.99 1:0.1:9.9 10:1:100].^n; xL=xR;
    %Recalculate pR and pL for new x-axis
    pL=[];pR=[];RHS_L=[];RHS_R=[];
    for i = 1:length(xL)
        for j = 1:length(xR)
            RHS_L(i,j) = B(1,1) + B(2,1)*xL(i) + B(3,1)*xR(j); 
            RHS_R(i,j) = B(1,2) + B(2,2)*xL(i) + B(3,2)*xR(j); 
            pR(i,j) = exp(RHS_R(i,j)) ./ (1+ exp(RHS_L(i,j)) + exp(RHS_R(i,j)));
        end
    end 
    xxR = xR.^(1/n);
    half = pR(1,1) + (pR(1,end)-pR(1,1))/2; %Half-detection value
    thresholds(iMuis).MTR = xxR(find(pR(1,:)>half,1)); %Visual Threshold Right
    modeldprimes(iMuis).MTR = norminv(pR(1,end))-norminv(pR(1,1));
end
%Tactile LEFT
for iMuis = 1:length(files)
    load(fullfile(files(iMuis).folder,files(iMuis).name)) 
    xR = [0:0.01:0.99 1:0.1:9.9 10:1:100].^n; xL=xR;
    %Recalculate pR and pL for new x-axis
    pL=[];pR=[];RHS_L=[];RHS_R=[];
    for i = 1:length(xL)
        for j = 1:length(xR)
            RHS_L(i,j) = B(1,1) + B(2,1)*xL(i) + B(3,1)*xR(j); 
            RHS_R(i,j) = B(1,2) + B(2,2)*xL(i) + B(3,2)*xR(j); 
            pL(i,j) = exp(RHS_L(i,j)) ./ (1+ exp(RHS_L(i,j)) + exp(RHS_R(i,j)));
        end
    end 
    xxL = xL.^(1/n);
    half = pL(1,1) + (pL(end,1)-pL(1,1))/2; %Half-detection value
    thresholds(iMuis).MTL = xxL(find(pL(:,1)>half,1)); %Visual Threshold Left
    modeldprimes(iMuis).MTL = norminv(pL(end,1))-norminv(pL(1,1));    
end

%% Stats
UV = [[thresholds.UVL] [thresholds.UVR]];
MV = [[thresholds.MVL] [thresholds.MVR]];
ranksum(UV,MV) %p=0.7018
%Try lme even if 2 points per mouse?
Threshold = [UV'; MV'];
MouseID = repmat([1 2 3 4]',4,1);
Cohort = [repmat({'UST'},8,1) ; repmat({'MST'},8,1)];
Side = [repmat({'L'},4,1) ; repmat({'R'},4,1) ; repmat({'L'},4,1) ; repmat({'R'},4,1)];
tbl = table(Threshold,MouseID,Cohort,Side);
tbl.MouseID = categorical(tbl.MouseID);
tbl.Cohort = categorical(tbl.Cohort);
tbl.Side = categorical(tbl.Side);
% lme = fitlme(tbl,'Threshold ~ Cohort + (1|MouseID)');
lme = fitlme(tbl,'Threshold ~ Cohort + Side');
stats = dataset2table(anova(lme,'DFMethod','Satterthwaite'));
p = stats{2,5};

% Plot
figure(); hold on
scatter(1*ones(size(UV)),UV, 'filled')
errorbar(2,nanmedian(UV),iqr(UV)/2,'o')
scatter(3*ones(size(MV)),MV, 'filled')
errorbar(4,nanmedian(MV),iqr(MV)/2,'o')
title('threshold')
xlim([0 5])
set(gca,'TickDir','out')
set(gca, 'YScale', 'log')
ylim([0.1 100])