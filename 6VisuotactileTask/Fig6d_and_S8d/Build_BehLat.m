%This code build the structure BehLat containing dprimes under different
%conditions

%% Get dprimes for Visual Contralateral
%% %%%%%%%Contralateral UST
folder = 'C:\Scratch\Beh\Data\Second Bump Unisensory - selection Jean 4';
files = dir(strcat(folder,'\*.mat'));
cc= 1 ;
BehLat = struct;

for iSess=1:length(files)
clearvars -except iSess files fi* ALL* cc BehLat
%Load session trialData (no need for neuronal activity)
load(fullfile(files(iSess).folder,files(iSess).name))

%Reparse and pre-process
pie_beh_reparse
trialData.stimDeflection(isnan(trialData.stimDeflection))=0;
trialData.stimContrast(isnan(trialData.stimContrast))=0;
if isfield(trialData,'stimVISSide')
    trialData.stimDeflection(cell2mat(trialData.stimTACSide)=='L') =  trialData.stimDeflection(cell2mat(trialData.stimTACSide)=='L')*-1;
    trialData.stimContrast(cell2mat(trialData.stimVISSide)=='L')   =  trialData.stimContrast(cell2mat(trialData.stimVISSide)=='L')*-1;    
else
    trialData.stimDeflection(trialData.leftCorrect==1) =  trialData.stimDeflection(trialData.leftCorrect==1)*-1;
    trialData.stimContrast(trialData.leftCorrect==1)   =  trialData.stimContrast(trialData.leftCorrect==1)*-1;
end
%Compute variables per condition
conds = 1:31; %unique(trialData.conditionnumber);
condsnummax = max(conds);
for k = 1:condsnummax
      TotalTrials(k) = nansum(trialData.conditionnumber==conds(k));  
      HitRight(k) = nansum(trialData.firstRespLFR=='R' & trialData.conditionnumber==conds(k)) / TotalTrials(k);
      HitLeft(k) = nansum(trialData.firstRespLFR=='L' & trialData.conditionnumber==conds(k)) / TotalTrials(k);
      NoResp(k)  = nansum(isnan(trialData.firstRespLFR) & trialData.conditionnumber==conds(k)) / TotalTrials(k);
      Contrasts(k) = mean(trialData.stimContrast(trialData.conditionnumber==conds(k)));
      Deflections(k) = mean(trialData.stimDeflection(trialData.conditionnumber==conds(k)));
      %CIM
      CorrectCh(k) = nansum(trialData.correctResponse==1 & trialData.conditionnumber==conds(k)) / TotalTrials(k);
      IncorrectCh(k) = nansum(trialData.firstIncorrect==1 & trialData.conditionnumber==conds(k)) / TotalTrials(k);
      NoRespCh(k) = nansum(trialData.noResponse==1 & trialData.conditionnumber==conds(k)) / TotalTrials(k);
end
%For Right-trials
righttrials = [1 2 3 4 8 9 12 13 16 17 18 19 22 23 26 27 28 29];
%Normalize HitRight !! ERASE PLS
% ginc = [repmat([6],1,15) repmat([21],1,21-16+1) 6 6 6 6 repmat([31],1,31-26+1)];
% linc = [repmat([3],1,15) repmat([18],1,21-16+1) 3 3 3 3 repmat([28],1,31-26+1)];
% for k = condsnummax
% HitRight(k) = (HitRight(k)-HitRight(ginc(k)))/(HitRight(linc(k)) - HitRight(ginc(k)));
% end
%Change %s so that they dont equal 0 or 100 and ruin the dprime. Add [1 0]
%to break
HitRight(HitRight==0) = 1./(TotalTrials(HitRight==0)+2);
HitRight(HitRight==1) = 1-1./(TotalTrials(HitRight==1)+2); 
HitLeft(HitLeft==0) = 1./(TotalTrials(HitLeft==0)+2);
HitLeft(HitLeft==1) = 1-1./(TotalTrials(HitLeft==1)+2);
%Compute dprimes per condition. 
DPR = [];
farincond = [ repmat([1],1,15) repmat([16],1,21-16+1) 1 1 1 1 repmat([26],1,31-26+1) ];

%dprime
for k = righttrials
    DPR(k) = norminv(HitRight(k)) - norminv(HitRight(farincond(k))) ;
    ALL_DPR(iSess,k) = DPR(k);
end

%Compute reduction in dprime per pairs of laser ON/OFF
DRED = [];
DRED(1) = ( DPR(2)-DPR(17) )/ DPR(2); %LOW+FULL
DRED(2) = ( DPR(3)-DPR(18) )/ DPR(3); %HIGH+FULL
DRED(3) = ( DPR(2)-DPR(27) )/ DPR(2); %LOW+DELAYED
DRED(4) = ( DPR(3)-DPR(28) )/ DPR(3); %HIGH+DELAYED

%Compute average RT
RELRT_DRED = [];
usedconds = unique(trialData.conditionnumber);
RT = cell(max(usedconds),2);
for k=1:length(trialData.stimStart)
    iCond = trialData.conditionnumber(k);
    LTIW = trialData.lickTime{k}(trialData.lickTime{k}>=trialData.respwinStart(k) & trialData.lickTime{k}<=(trialData.respwinEnd(k)+0.002));
    if ~isempty(LTIW) && iCond~=0 
        if trialData.correctResponse(k)
        RT{iCond,1} = [ RT{iCond,1} (LTIW(1) - trialData.stimStart(k)) ];
        else %elseif trialData.firstIncorrect(k)
        RT{iCond,2} = [ RT{iCond,2} (LTIW(1) - trialData.stimStart(k)) ];    
        end
    end
end
clear ans iCond k LTIW usedconds
RELRT_DRED(1) = 0 - mean(RT{2,1}); %full-low - use the control condition 
RELRT_DRED(2) = 0 - mean(RT{3,1}); %full-high
RELRT_DRED(3) = 0.240 - mean(RT{2,1}); %delayed-low
RELRT_DRED(4) = 0.240 - mean(RT{3,1}); %delayed-high

%Do same (per condition) in struct
%Cond1
BehLat(cc).Task             = 'tUni';       %tUni/tMulti                
BehLat(cc).Saliency         = 'Low';        %Low/High
BehLat(cc).Silencing        = 'Early';      %Off/Early/Late
BehLat(cc).Modality         = 'Visual';      
BehLat(cc).Side             = 'Contra';  
BehLat(cc).LaserOnset       = 0;            %Nan/0/0.240 (in s)
BehLat(cc).Dprime           = DPR(17);
BehLat(cc).DprimeCtrl       = DPR(2);
BehLat(cc).DprimeReduction  = ( DPR(2)-DPR(17) )/ DPR(2);
BehLat(cc).RTmean           = mean(RT{17,1});
BehLat(cc).RTctrl           = mean(RT{2,1});
BehLat(cc).RTctrlmed        = median(RT{2,1});
BehLat(cc).SessionName      = files(iSess).name;
%
BehLat(cc+1).Task             = 'tUni';       %tUni/tMulti                
BehLat(cc+1).Saliency         = 'High';        %Low/High
BehLat(cc+1).Silencing        = 'Early';      %Off/Early/Late
BehLat(cc+1).Modality         = 'Visual';      
BehLat(cc+1).Side             = 'Contra';  
BehLat(cc+1).LaserOnset       = 0;            %Nan/0/0.240 (in s)
BehLat(cc+1).Dprime           = DPR(18);
BehLat(cc+1).DprimeCtrl       = DPR(3);
BehLat(cc+1).DprimeReduction  = ( DPR(3)-DPR(18) )/ DPR(3);
BehLat(cc+1).RTmean           = mean(RT{18,1});
BehLat(cc+1).RTctrl           = mean(RT{3,1});
BehLat(cc+1).RTctrlmed        = median(RT{3,1});
BehLat(cc+1).SessionName      = files(iSess).name;
%
BehLat(cc+2).Task             = 'tUni';       %tUni/tMulti                
BehLat(cc+2).Saliency         = 'Low';        %Low/High
BehLat(cc+2).Silencing        = 'Late';      %Off/Early/Late
BehLat(cc+2).Modality         = 'Visual';      
BehLat(cc+2).Side             = 'Contra';  
BehLat(cc+2).LaserOnset       = 0.240;            %Nan/0/0.240 (in s)
BehLat(cc+2).Dprime           = DPR(27);
BehLat(cc+2).DprimeCtrl       = DPR(2);
BehLat(cc+2).DprimeReduction  = ( DPR(2)-DPR(27) )/ DPR(2);
BehLat(cc+2).RTmean           = mean(RT{27,1});
BehLat(cc+2).RTctrl           = mean(RT{2,1});
BehLat(cc+2).RTctrlmed        = median(RT{2,1});
BehLat(cc+2).SessionName      = files(iSess).name;
%
BehLat(cc+3).Task             = 'tUni';       %tUni/tMulti                
BehLat(cc+3).Saliency         = 'High';        %Low/High
BehLat(cc+3).Silencing        = 'Late';      %Off/Early/Late
BehLat(cc+3).Modality         = 'Visual';      
BehLat(cc+3).Side             = 'Contra';  
BehLat(cc+3).LaserOnset       = 0.240;            %Nan/0/0.240 (in s)
BehLat(cc+3).Dprime           = DPR(28);
BehLat(cc+3).DprimeCtrl       = DPR(3);
BehLat(cc+3).DprimeReduction  = ( DPR(3)-DPR(28) )/ DPR(3);
BehLat(cc+3).RTmean           = mean(RT{28,1});
BehLat(cc+3).RTctrl           = mean(RT{3,1});
BehLat(cc+3).RTctrlmed        = median(RT{3,1});
BehLat(cc+3).SessionName      = files(iSess).name;
%
cc = cc + 4;
end
clearvars -except BehLat
%% %%%%%%%Contralateral MST
folder_full = 'C:\Scratch\Beh\Data\V1 second bump silencing\V1 full stim silencing';
folder_delayed = 'C:\Scratch\Beh\Data\V1 second bump silencing\V1 delayed silencing';
files_full = dir(strcat(folder_full,'\*.mat'));
files_delayed = dir(strcat(folder_delayed,'\*.mat'));
cc= length(BehLat)+1 ;
% FOR FULL
for iSess=1:length(files_full)
clearvars -except iSess files* fi* ALL* cc BehLat
%Load session trialData (no need for neuronal activity)
load(fullfile(files_full(iSess).folder,files_full(iSess).name))

%Reparse and pre-process
pie_beh_reparse
trialData.stimDeflection(isnan(trialData.stimDeflection))=0;
trialData.stimContrast(isnan(trialData.stimContrast))=0;
if isfield(trialData,'stimVISSide')
    trialData.stimDeflection(cell2mat(trialData.stimTACSide)=='L') =  trialData.stimDeflection(cell2mat(trialData.stimTACSide)=='L')*-1;
    trialData.stimContrast(cell2mat(trialData.stimVISSide)=='L')   =  trialData.stimContrast(cell2mat(trialData.stimVISSide)=='L')*-1;    
else
    trialData.stimDeflection(trialData.leftCorrect==1) =  trialData.stimDeflection(trialData.leftCorrect==1)*-1;
    trialData.stimContrast(trialData.leftCorrect==1)   =  trialData.stimContrast(trialData.leftCorrect==1)*-1;
end
%Compute variables per condition
conds = 1:42; %unique(trialData.conditionnumber);
condsnummax = max(conds);
usedconds = unique(trialData.conditionnumber);
for k = 1:condsnummax
      TotalTrials(k) = nansum(trialData.conditionnumber==conds(k));  
      HitRight(k) = nansum(trialData.firstRespLFR=='R' & trialData.conditionnumber==conds(k)) / TotalTrials(k);
      HitLeft(k) = nansum(trialData.firstRespLFR=='L' & trialData.conditionnumber==conds(k)) / TotalTrials(k);
      NoResp(k)  = nansum(isnan(trialData.firstRespLFR) & trialData.conditionnumber==conds(k)) / TotalTrials(k);
      Contrasts(k) = mean(trialData.stimContrast(trialData.conditionnumber==conds(k)));
      Deflections(k) = mean(trialData.stimDeflection(trialData.conditionnumber==conds(k)));
      %CIM
      CorrectCh(k) = nansum(trialData.correctResponse==1 & trialData.conditionnumber==conds(k)) / TotalTrials(k);
      IncorrectCh(k) = nansum(trialData.firstIncorrect==1 & trialData.conditionnumber==conds(k)) / TotalTrials(k);
      NoRespCh(k) = nansum(trialData.noResponse==1 & trialData.conditionnumber==conds(k)) / TotalTrials(k);
end
%Change %s so that they dont equal 0 or 100 and ruin the dprime
HitRight(HitRight==0) = 1./(TotalTrials(HitRight==0)+2);
HitRight(HitRight==1) = 1-1./(TotalTrials(HitRight==1)+2);
HitLeft(HitLeft==0) = 1./(TotalTrials(HitLeft==0)+2);
HitLeft(HitLeft==1) = 1-1./(TotalTrials(HitLeft==1)+2);
%Compute dprimes per condition. 
DPR = [];
righttrials = [1  2  3  4  8  9 12 13 16 17 18 19 20 21 29 30 31 32 37 38 39 40 41 42];
%farincond uses normal absolute condition number as index !
farincond   = [repmat(1,1,15) repmat(16,1,21-16+1) 1 1 1 1 16 16 16 1 1 1 1 1 1 1 1 16 16 16 16 16 16];
for k = righttrials(ismember(righttrials,usedconds))
    DPR(k) = norminv(HitRight(k)) - norminv(HitRight(farincond(k))) ;
end
%Compute reduction in dprime per pairs of laser ON/OFF
DRED = [];
DRED(1) = ( DPR(2)-DPR(17) )/ DPR(2); %LOW+FULL
DRED(2) = ( DPR(3)-DPR(18) )/ DPR(3); %HIGH+FULL
% DRED(3) = ( DPR(2)-DPR(27) )/ DPR(2); %LOW+FULL MULTI
% DRED(4) = ( DPR(3)-DPR(28) )/ DPR(3); %HIGH+FULL  MULTI

%Compute average RT
RELRT_DRED = [];
usedconds = unique(trialData.conditionnumber);
RT = cell(max(usedconds),2);
for k=1:length(trialData.stimStart)
    iCond = trialData.conditionnumber(k);
    LTIW = trialData.lickTime{k}(trialData.lickTime{k}>=trialData.respwinStart(k) & trialData.lickTime{k}<=(trialData.respwinEnd(k)+0.002));
    if ~isempty(LTIW) && iCond~=0 
        if trialData.correctResponse(k)
        RT{iCond,1} = [ RT{iCond,1} (LTIW(1) - trialData.stimStart(k)) ];
        else %elseif trialData.firstIncorrect(k)
        RT{iCond,2} = [ RT{iCond,2} (LTIW(1) - trialData.stimStart(k)) ];    
        end
    end
end
clear ans iCond k LTIW usedconds
RELRT_DRED(1) = 0       - mean(RT{2,1}); %full-low
RELRT_DRED(2) = 0       - mean(RT{3,1}); %full-high
%RELRT_DRED(3) = 0      - mean(RT{12,1}); %full-low MULTI
%RELRT_DRED(4) = 0      - mean(RT{13,1}); %full-high MULTI

%Fill BehLatencies struct per condition
BehLat(cc).Task             = 'tMulti';       %tUni/tMulti                
BehLat(cc).Saliency         = 'Low';        %Low/High
BehLat(cc).Silencing        = 'Early';      %Off/Early/Late
BehLat(cc).Modality         = 'Visual';      
BehLat(cc).Side             = 'Contra';  
BehLat(cc).LaserOnset       = 0;            %Nan/0/0.240 (in s)
BehLat(cc).Dprime           = DPR(17);
BehLat(cc).DprimeCtrl       = DPR(2);
BehLat(cc).DprimeReduction  = ( DPR(2)-DPR(17) )/ DPR(2);
BehLat(cc).RTmean           = mean(RT{17,1});
BehLat(cc).RTctrl           = mean(RT{2,1});
BehLat(cc).RTctrlmed        = median(RT{2,1});
BehLat(cc).SessionName      = files_full(iSess).name;
%
BehLat(cc+1).Task             = 'tMulti';       %tUni/tMulti                
BehLat(cc+1).Saliency         = 'High';        %Low/High
BehLat(cc+1).Silencing        = 'Early';      %Off/Early/Late
BehLat(cc+1).Modality         = 'Visual';      
BehLat(cc+1).Side             = 'Contra';  
BehLat(cc+1).LaserOnset       = 0;            %Nan/0/0.240 (in s)
BehLat(cc+1).Dprime           = DPR(18);
BehLat(cc+1).DprimeCtrl       = DPR(3);
BehLat(cc+1).DprimeReduction  = ( DPR(3)-DPR(18) )/ DPR(3);
BehLat(cc+1).RTmean           = mean(RT{18,1});
BehLat(cc+1).RTctrl           = mean(RT{3,1});
BehLat(cc+1).RTctrlmed        = median(RT{3,1});
BehLat(cc+1).SessionName      = files_full(iSess).name;
%
cc= cc+2;
end
% FOR DELAYED
for iSess=1:length(files_delayed)
clearvars -except iSess files* fi* ALL* cc BehLat
%Load session trialData (no need for neuronal activity)
load(fullfile(files_delayed(iSess).folder,files_delayed(iSess).name))

%Reparse and pre-process
pie_beh_reparse
trialData.stimDeflection(isnan(trialData.stimDeflection))=0;
trialData.stimContrast(isnan(trialData.stimContrast))=0;
if isfield(trialData,'stimVISSide')
    trialData.stimDeflection(cell2mat(trialData.stimTACSide)=='L') =  trialData.stimDeflection(cell2mat(trialData.stimTACSide)=='L')*-1;
    trialData.stimContrast(cell2mat(trialData.stimVISSide)=='L')   =  trialData.stimContrast(cell2mat(trialData.stimVISSide)=='L')*-1;    
else
    trialData.stimDeflection(trialData.leftCorrect==1) =  trialData.stimDeflection(trialData.leftCorrect==1)*-1;
    trialData.stimContrast(trialData.leftCorrect==1)   =  trialData.stimContrast(trialData.leftCorrect==1)*-1;
end
%Compute variables per condition
conds = 1:42; %unique(trialData.conditionnumber);
condsnummax = max(conds);
usedconds = unique(trialData.conditionnumber);
for k = 1:condsnummax
      TotalTrials(k) = nansum(trialData.conditionnumber==conds(k));  
      HitRight(k) = nansum(trialData.firstRespLFR=='R' & trialData.conditionnumber==conds(k)) / TotalTrials(k);
      HitLeft(k) = nansum(trialData.firstRespLFR=='L' & trialData.conditionnumber==conds(k)) / TotalTrials(k);
      NoResp(k)  = nansum(isnan(trialData.firstRespLFR) & trialData.conditionnumber==conds(k)) / TotalTrials(k);
      Contrasts(k) = mean(trialData.stimContrast(trialData.conditionnumber==conds(k)));
      Deflections(k) = mean(trialData.stimDeflection(trialData.conditionnumber==conds(k)));
      %CIM
      CorrectCh(k) = nansum(trialData.correctResponse==1 & trialData.conditionnumber==conds(k)) / TotalTrials(k);
      IncorrectCh(k) = nansum(trialData.firstIncorrect==1 & trialData.conditionnumber==conds(k)) / TotalTrials(k);
      NoRespCh(k) = nansum(trialData.noResponse==1 & trialData.conditionnumber==conds(k)) / TotalTrials(k);
end
%Change %s so that they dont equal 0 or 100 and ruin the dprime
HitRight(HitRight==0) = 1./(TotalTrials(HitRight==0)+2);
HitRight(HitRight==1) = 1-1./(TotalTrials(HitRight==1)+2);
HitLeft(HitLeft==0) = 1./(TotalTrials(HitLeft==0)+2);
HitLeft(HitLeft==1) = 1-1./(TotalTrials(HitLeft==1)+2);
%Compute dprimes per condition. 
DPR = [];
righttrials = [1  2  3  4  8  9 12 13 16 17 18 19 20 21 29 30 31 32 37 38 39 40 41 42];
%farincond uses normal absolute condition number as index !
farincond   = [repmat(1,1,15) repmat(16,1,21-16+1) 1 1 1 1 16 16 16 1 1 1 1 1 1 1 1 16 16 16 16 16 16];
for k = righttrials(ismember(righttrials,usedconds))
    DPR(k) = norminv(HitRight(k)) - norminv(HitRight(farincond(k))) ;
end
%Compute reduction in dprime per pairs of laser ON/OFF
DRED = [];
DRED(1) = ( DPR(2)-DPR(17) )/ DPR(2); %LOW+FULL
DRED(2) = ( DPR(3)-DPR(18) )/ DPR(3); %HIGH+FULL
% DRED(3) = ( DPR(2)-DPR(27) )/ DPR(2); %LOW+FULL MULTI
% DRED(4) = ( DPR(3)-DPR(28) )/ DPR(3); %HIGH+FULL  MULTI

%Compute average RT
RELRT_DRED = [];
usedconds = unique(trialData.conditionnumber);
RT = cell(max(usedconds),2);
for k=1:length(trialData.stimStart)
    iCond = trialData.conditionnumber(k);
    LTIW = trialData.lickTime{k}(trialData.lickTime{k}>=trialData.respwinStart(k) & trialData.lickTime{k}<=(trialData.respwinEnd(k)+0.002));
    if ~isempty(LTIW) && iCond~=0 
        if trialData.correctResponse(k)
        RT{iCond,1} = [ RT{iCond,1} (LTIW(1) - trialData.stimStart(k)) ];
        else %elseif trialData.firstIncorrect(k)
        RT{iCond,2} = [ RT{iCond,2} (LTIW(1) - trialData.stimStart(k)) ];    
        end
    end
end
clear ans iCond k LTIW usedconds
RELRT_DRED(1) = 0.240       - mean(RT{2,1}); %delayed-low
RELRT_DRED(2) = 0.240       - mean(RT{3,1}); %delayed-high
%RELRT_DRED(3) = 0      - mean(RT{12,1}); %delayed-low MULTI
%RELRT_DRED(4) = 0      - mean(RT{13,1}); %delayed-high MULTI

%Fill BehLatencies struct per condition
BehLat(cc).Task             = 'tMulti';       %tUni/tMulti                
BehLat(cc).Saliency         = 'Low';        %Low/High
BehLat(cc).Silencing        = 'Late';      %Off/Early/Late
BehLat(cc).Modality         = 'Visual';      
BehLat(cc).Side             = 'Contra';  
BehLat(cc).LaserOnset       = 0.240;            %Nan/0/0.240 (in s)
BehLat(cc).Dprime           = DPR(17);
BehLat(cc).DprimeCtrl       = DPR(2);
BehLat(cc).DprimeReduction  = ( DPR(2)-DPR(17) )/ DPR(2);
BehLat(cc).RTmean           = mean(RT{17,1});
BehLat(cc).RTctrl           = mean(RT{2,1});
BehLat(cc).RTctrlmed        = median(RT{2,1});
BehLat(cc).SessionName      = files_delayed(iSess).name;
%
BehLat(cc+1).Task             = 'tMulti';       %tUni/tMulti                
BehLat(cc+1).Saliency         = 'High';        %Low/High
BehLat(cc+1).Silencing        = 'Late';      %Off/Early/Late
BehLat(cc+1).Modality         = 'Visual';      
BehLat(cc+1).Side             = 'Contra';  
BehLat(cc+1).LaserOnset       = 0.240;            %Nan/0/0.240 (in s)
BehLat(cc+1).Dprime           = DPR(18);
BehLat(cc+1).DprimeCtrl       = DPR(3);
BehLat(cc+1).DprimeReduction  = ( DPR(3)-DPR(18) )/ DPR(3);
BehLat(cc+1).RTmean           = mean(RT{18,1});
BehLat(cc+1).RTctrl           = mean(RT{3,1});
BehLat(cc+1).RTctrlmed        = median(RT{3,1});
BehLat(cc+1).SessionName      = files_delayed(iSess).name;
%
cc= cc+2;
end
clearvars -except BehLat 

%% Get dprimes for Visual Ipsilateral 
%% %%%%%%%Ipsilateral UST
folder = 'C:\Scratch\Beh\Data\Second Bump Unisensory - selection Jean 4';
files = dir(strcat(folder,'\*.mat'));
cc= length(BehLat)+1 ;
for iSess=1:length(files)
clearvars -except iSess files fi* ALL* cc BehLat
%Load session trialData (no need for neuronal activity)
load(fullfile(files(iSess).folder,files(iSess).name))

%Reparse and pre-process
pie_beh_reparse
trialData.stimDeflection(isnan(trialData.stimDeflection))=0;
trialData.stimContrast(isnan(trialData.stimContrast))=0;
if isfield(trialData,'stimVISSide')
    trialData.stimDeflection(cell2mat(trialData.stimTACSide)=='L') =  trialData.stimDeflection(cell2mat(trialData.stimTACSide)=='L')*-1;
    trialData.stimContrast(cell2mat(trialData.stimVISSide)=='L')   =  trialData.stimContrast(cell2mat(trialData.stimVISSide)=='L')*-1;    
else
    trialData.stimDeflection(trialData.leftCorrect==1) =  trialData.stimDeflection(trialData.leftCorrect==1)*-1;
    trialData.stimContrast(trialData.leftCorrect==1)   =  trialData.stimContrast(trialData.leftCorrect==1)*-1;
end
%Compute variables per condition
conds = 1:31; %unique(trialData.conditionnumber);
condsnummax = max(conds);
usedconds = unique(trialData.conditionnumber);
for k = 1:condsnummax
      TotalTrials(k) = nansum(trialData.conditionnumber==conds(k));  
      HitRight(k) = nansum(trialData.firstRespLFR=='R' & trialData.conditionnumber==conds(k)) / TotalTrials(k);
      HitLeft(k) = nansum(trialData.firstRespLFR=='L' & trialData.conditionnumber==conds(k)) / TotalTrials(k);
      NoResp(k)  = nansum(isnan(trialData.firstRespLFR) & trialData.conditionnumber==conds(k)) / TotalTrials(k);
      Contrasts(k) = mean(trialData.stimContrast(trialData.conditionnumber==conds(k)));
      Deflections(k) = mean(trialData.stimDeflection(trialData.conditionnumber==conds(k)));
end
%For Right-trials
lefttrials = [1 5 6 7 10 11 14 15 16 20 21 24 25 26 30 31];

HitRight(HitRight==0) = 1./(TotalTrials(HitRight==0)+2);
HitRight(HitRight==1) = 1-1./(TotalTrials(HitRight==1)+2); 
HitLeft(HitLeft==0) = 1./(TotalTrials(HitLeft==0)+2);
HitLeft(HitLeft==1) = 1-1./(TotalTrials(HitLeft==1)+2);
%Compute dprimes per condition. 
DPR = [];
%farincond must match all conds not just selected ones
farincond = [ repmat([1],1,15) repmat([16],1,21-16+1) 1 1 1 1 repmat([26],1,31-26+1) ];

%dprime
for k = lefttrials(ismember(lefttrials,usedconds))
    DPR(k) = norminv(HitLeft(k)) - norminv(HitLeft(farincond(k))) ;
    ALL_DPR(iSess,k) = DPR(k);
end

%Compute reduction in dprime per pairs of laser ON/OFF
DRED = [];
DRED(1) = ( DPR(5)-DPR(20) )/ DPR(5); %LOW+FULL
DRED(2) = ( DPR(6)-DPR(21) )/ DPR(21); %HIGH+FULL
DRED(3) = ( DPR(5)-DPR(30) )/ DPR(5); %LOW+DELAYED
DRED(4) = ( DPR(6)-DPR(31) )/ DPR(6); %HIGH+DELAYED

%Compute average RT
RELRT_DRED = [];
usedconds = unique(trialData.conditionnumber);
RT = cell(max(usedconds),2);
for k=1:length(trialData.stimStart)
    iCond = trialData.conditionnumber(k);
    LTIW = trialData.lickTime{k}(trialData.lickTime{k}>=trialData.respwinStart(k) & trialData.lickTime{k}<=(trialData.respwinEnd(k)+0.002));
    if ~isempty(LTIW) && iCond~=0 
        if trialData.correctResponse(k)
        RT{iCond,1} = [ RT{iCond,1} (LTIW(1) - trialData.stimStart(k)) ];
        else %elseif trialData.firstIncorrect(k)
        RT{iCond,2} = [ RT{iCond,2} (LTIW(1) - trialData.stimStart(k)) ];    
        end
    end
end
clear ans iCond k LTIW usedconds
RELRT_DRED(1) = 0 - mean(RT{5,1}); %full-low - use the control condition 
RELRT_DRED(2) = 0 - mean(RT{6,1}); %full-high
RELRT_DRED(3) = 0.240 - mean(RT{5,1}); %delayed-low
RELRT_DRED(4) = 0.240 - mean(RT{6,1}); %delayed-high

%Do same (per condition) in struct
%Cond1
BehLat(cc).Task             = 'tUni';       %tUni/tMulti                
BehLat(cc).Saliency         = 'Low';        %Low/High
BehLat(cc).Silencing        = 'Early';      %Off/Early/Late
BehLat(cc).Modality         = 'Visual';
BehLat(cc).Side             = 'Ipsi';
BehLat(cc).LaserOnset       = 0;            %Nan/0/0.240 (in s)
BehLat(cc).Dprime           = DPR(20);
BehLat(cc).DprimeCtrl       = DPR(5);
BehLat(cc).DprimeReduction  = ( DPR(5)-DPR(20) )/ DPR(5);
BehLat(cc).RTmean           = mean(RT{20,1});
BehLat(cc).RTctrl           = mean(RT{5,1});
BehLat(cc).RTctrlmed        = median(RT{5,1});
BehLat(cc).SessionName      = files(iSess).name;
%
BehLat(cc+1).Task             = 'tUni';       %tUni/tMulti                
BehLat(cc+1).Saliency         = 'High';        %Low/High
BehLat(cc+1).Silencing        = 'Early';      %Off/Early/Late
BehLat(cc+1).Modality         = 'Visual';
BehLat(cc+1).Side             = 'Ipsi';
BehLat(cc+1).LaserOnset       = 0;            %Nan/0/0.240 (in s)
BehLat(cc+1).Dprime           = DPR(21);
BehLat(cc+1).DprimeCtrl       = DPR(6);
BehLat(cc+1).DprimeReduction  = ( DPR(6)-DPR(21) )/ DPR(6);
BehLat(cc+1).RTmean           = mean(RT{21,1});
BehLat(cc+1).RTctrl           = mean(RT{6,1});
BehLat(cc+1).RTctrlmed        = median(RT{6,1});
BehLat(cc+1).SessionName      = files(iSess).name;
%
BehLat(cc+2).Task             = 'tUni';       %tUni/tMulti                
BehLat(cc+2).Saliency         = 'Low';        %Low/High
BehLat(cc+2).Silencing        = 'Late';      %Off/Early/Late
BehLat(cc+2).Modality         = 'Visual';
BehLat(cc+2).Side             = 'Ipsi';
BehLat(cc+2).LaserOnset       = 0.240;            %Nan/0/0.240 (in s)
BehLat(cc+2).Dprime           = DPR(30);
BehLat(cc+2).DprimeCtrl       = DPR(5);
BehLat(cc+2).DprimeReduction  = ( DPR(5)-DPR(30) )/ DPR(5);
BehLat(cc+2).RTmean           = mean(RT{30,1});
BehLat(cc+2).RTctrl           = mean(RT{5,1});
BehLat(cc+2).RTctrlmed        = median(RT{5,1});
BehLat(cc+2).SessionName      = files(iSess).name;
%
BehLat(cc+3).Task             = 'tUni';       %tUni/tMulti                
BehLat(cc+3).Saliency         = 'High';        %Low/High
BehLat(cc+3).Silencing        = 'Late';      %Off/Early/Late
BehLat(cc+3).Modality         = 'Visual';
BehLat(cc+3).Side             = 'Ipsi';
BehLat(cc+3).LaserOnset       = 0.240;            %Nan/0/0.240 (in s)
BehLat(cc+3).Dprime           = DPR(31);
BehLat(cc+3).DprimeCtrl       = DPR(6);
BehLat(cc+3).DprimeReduction  = ( DPR(6)-DPR(31) )/ DPR(6);
BehLat(cc+3).RTmean           = mean(RT{31,1});
BehLat(cc+3).RTctrl           = mean(RT{6,1});
BehLat(cc+3).RTctrlmed        = median(RT{6,1});
BehLat(cc+3).SessionName      = files(iSess).name;
%
cc = cc + 4;
end
clearvars -except BehLat

%% %%%%%%%Ipsilateral MST
folder_full = 'C:\Scratch\Beh\Data\V1 second bump silencing\V1 full stim silencing';
folder_delayed = 'C:\Scratch\Beh\Data\V1 second bump silencing\V1 delayed silencing';
files_full = dir(strcat(folder_full,'\*.mat'));
files_delayed = dir(strcat(folder_delayed,'\*.mat'));
cc= length(BehLat)+1 ;
% FOR FULL
for iSess=1:length(files_full)
clearvars -except iSess files* fi* ALL* cc BehLat
%Load session trialData (no need for neuronal activity)
load(fullfile(files_full(iSess).folder,files_full(iSess).name))

%Reparse and pre-process
pie_beh_reparse
trialData.stimDeflection(isnan(trialData.stimDeflection))=0;
trialData.stimContrast(isnan(trialData.stimContrast))=0;
if isfield(trialData,'stimVISSide')
    trialData.stimDeflection(cell2mat(trialData.stimTACSide)=='L') =  trialData.stimDeflection(cell2mat(trialData.stimTACSide)=='L')*-1;
    trialData.stimContrast(cell2mat(trialData.stimVISSide)=='L')   =  trialData.stimContrast(cell2mat(trialData.stimVISSide)=='L')*-1;    
else
    trialData.stimDeflection(trialData.leftCorrect==1) =  trialData.stimDeflection(trialData.leftCorrect==1)*-1;
    trialData.stimContrast(trialData.leftCorrect==1)   =  trialData.stimContrast(trialData.leftCorrect==1)*-1;
end

%Determine if OptoLR or not (if not then conditions 40 and 42 are multi
%right and not visual left
formatt = 'undefined';
if sum(trialData.conditionnumber==40)>0
   if mean(trialData.leftCorrect(trialData.conditionnumber==40))==1
       formatt= 'Rec2LR';
   else
       formatt= 'Rec2';
   end
elseif sum(trialData.conditionnumber==42)>0
   if mean(trialData.leftCorrect(trialData.conditionnumber==42))==1
       formatt= 'Rec2LR';
   else 
       formatt= 'Rec2';
   end
else %cant determine/doesn't matter
end

%Compute variables per condition
conds = 1:42; %unique(trialData.conditionnumber);
condsnummax = max(conds);
usedconds = unique(trialData.conditionnumber);
for k = 1:condsnummax
      TotalTrials(k) = nansum(trialData.conditionnumber==conds(k));  
      HitRight(k) = nansum(trialData.firstRespLFR=='R' & trialData.conditionnumber==conds(k)) / TotalTrials(k);
      HitLeft(k) = nansum(trialData.firstRespLFR=='L' & trialData.conditionnumber==conds(k)) / TotalTrials(k);
      NoResp(k)  = nansum(isnan(trialData.firstRespLFR) & trialData.conditionnumber==conds(k)) / TotalTrials(k);
      Contrasts(k) = mean(trialData.stimContrast(trialData.conditionnumber==conds(k)));
      Deflections(k) = mean(trialData.stimDeflection(trialData.conditionnumber==conds(k)));
end
%Change %s so that they dont equal 0 or 100 and ruin the dprime
HitRight(HitRight==0) = 1./(TotalTrials(HitRight==0)+2);
HitRight(HitRight==1) = 1-1./(TotalTrials(HitRight==1)+2);
HitLeft(HitLeft==0) = 1./(TotalTrials(HitLeft==0)+2);
HitLeft(HitLeft==1) = 1-1./(TotalTrials(HitLeft==1)+2);
%Compute dprimes per condition. 
DPR = [];
lefttrials = [1 5 6 7 10 11 14 15 40 42];
%farincond uses normal absolute condition number as index !
farincond   = [repmat(1,1,15) repmat(16,1,21-16+1) 1 1 1 1 16 16 16 1 1 1 1 1 1 1 1 16 16 16 16 16 16];
for k = lefttrials(ismember(lefttrials,usedconds))
    DPR(k) = norminv(HitLeft(k)) - norminv(HitLeft(farincond(k))) ;
end
%Compute reduction in dprime per pairs of laser ON/OFF
DRED = [];
if strcmp(formatt,'Rec2LR')
DRED(1) = ( DPR(5)-DPR(40) )/ DPR(5); %LOW+FULL
else
DRED(1) = nan; %LOW+FULL 
end
DRED(2) = nan; %Never tested
% DRED(3) = ( DPR(2)-DPR(27) )/ DPR(2); %LOW+FULL MULTI
% DRED(4) = ( DPR(3)-DPR(28) )/ DPR(3); %HIGH+FULL  MULTI

%Compute average RT
RELRT_DRED = [];
usedconds = unique(trialData.conditionnumber);
RT = cell(max(usedconds),2);
for k=1:length(trialData.stimStart)
    iCond = trialData.conditionnumber(k);
    LTIW = trialData.lickTime{k}(trialData.lickTime{k}>=trialData.respwinStart(k) & trialData.lickTime{k}<=(trialData.respwinEnd(k)+0.002));
    if ~isempty(LTIW) && iCond~=0 
        if trialData.correctResponse(k)
        RT{iCond,1} = [ RT{iCond,1} (LTIW(1) - trialData.stimStart(k)) ];
        else %elseif trialData.firstIncorrect(k)
        RT{iCond,2} = [ RT{iCond,2} (LTIW(1) - trialData.stimStart(k)) ];    
        end
    end
end
clear ans iCond k LTIW usedconds
RELRT_DRED(1) = 0       - mean(RT{5,1}); %full-low
RELRT_DRED(2) = 0       - mean(RT{6,1}); %full-high
%RELRT_DRED(3) = 0      - mean(RT{12,1}); %full-low MULTI
%RELRT_DRED(4) = 0      - mean(RT{13,1}); %full-high MULTI

%Fill BehLatencies struct per condition
BehLat(cc).Task             = 'tMulti';       %tUni/tMulti                
BehLat(cc).Saliency         = 'Low';        %Low/High
BehLat(cc).Silencing        = 'Early';      %Off/Early/Late
BehLat(cc).Modality         = 'Visual';
BehLat(cc).Side             = 'Ipsi';      
BehLat(cc).LaserOnset       = 0;            %Nan/0/0.240 (in s)
BehLat(cc).DprimeCtrl       = DPR(5);
if strcmp(formatt,'Rec2LR')
BehLat(cc).Dprime           = DPR(40);
BehLat(cc).DprimeReduction  = ( DPR(5)-DPR(40) )/ DPR(5);
BehLat(cc).RTmean           = mean(RT{40,1});
else
BehLat(cc).Dprime           = nan;
BehLat(cc).DprimeReduction  = nan;
BehLat(cc).RTmean           = nan;
end
BehLat(cc).RTctrl           = mean(RT{5,1});
BehLat(cc).RTctrlmed        = median(RT{5,1});
BehLat(cc).SessionName      = files_full(iSess).name;
%
BehLat(cc+1).Task             = 'tMulti';       %tUni/tMulti                
BehLat(cc+1).Saliency         = 'High';        %Low/High
BehLat(cc+1).Silencing        = 'Early';      %Off/Early/Late
BehLat(cc+1).Modality         = 'Visual';
BehLat(cc+1).Side             = 'Ipsi';      
BehLat(cc+1).LaserOnset       = 0;            %Nan/0/0.240 (in s)
BehLat(cc+1).Dprime           = nan;
BehLat(cc+1).DprimeCtrl       = DPR(6);
BehLat(cc+1).DprimeReduction  = nan;
BehLat(cc+1).RTmean           = nan;
BehLat(cc+1).RTctrl           = mean(RT{6,1});
BehLat(cc+1).RTctrlmed        = median(RT{6,1});
BehLat(cc+1).SessionName      = files_full(iSess).name;
%
cc= cc+2;
end
% FOR DELAYED
for iSess=1:length(files_delayed)
clearvars -except iSess files* fi* ALL* cc BehLat
%Load session trialData (no need for neuronal activity)
load(fullfile(files_delayed(iSess).folder,files_delayed(iSess).name))

%Reparse and pre-process
pie_beh_reparse
trialData.stimDeflection(isnan(trialData.stimDeflection))=0;
trialData.stimContrast(isnan(trialData.stimContrast))=0;
if isfield(trialData,'stimVISSide')
    trialData.stimDeflection(cell2mat(trialData.stimTACSide)=='L') =  trialData.stimDeflection(cell2mat(trialData.stimTACSide)=='L')*-1;
    trialData.stimContrast(cell2mat(trialData.stimVISSide)=='L')   =  trialData.stimContrast(cell2mat(trialData.stimVISSide)=='L')*-1;    
else
    trialData.stimDeflection(trialData.leftCorrect==1) =  trialData.stimDeflection(trialData.leftCorrect==1)*-1;
    trialData.stimContrast(trialData.leftCorrect==1)   =  trialData.stimContrast(trialData.leftCorrect==1)*-1;
end
%Determine if OptoLR or not (if not then conditions 40 and 42 are multi
%right and not visual left
formatt = 'undefined';
if sum(trialData.conditionnumber==40)>0
   if mean(trialData.leftCorrect(trialData.conditionnumber==40))==1
       formatt= 'Rec2LR';
   else
       formatt= 'Rec2';
   end
elseif sum(trialData.conditionnumber==42)>0
   if mean(trialData.leftCorrect(trialData.conditionnumber==42))==1
       formatt= 'Rec2LR';
   else 
       formatt= 'Rec2';
   end
else %cant determine/doesn't matter
end
%Compute variables per condition
conds = 1:42; %unique(trialData.conditionnumber);
condsnummax = max(conds);
usedconds = unique(trialData.conditionnumber);
for k = 1:condsnummax
      TotalTrials(k) = nansum(trialData.conditionnumber==conds(k));  
      HitRight(k) = nansum(trialData.firstRespLFR=='R' & trialData.conditionnumber==conds(k)) / TotalTrials(k);
      HitLeft(k) = nansum(trialData.firstRespLFR=='L' & trialData.conditionnumber==conds(k)) / TotalTrials(k);
      NoResp(k)  = nansum(isnan(trialData.firstRespLFR) & trialData.conditionnumber==conds(k)) / TotalTrials(k);
      Contrasts(k) = mean(trialData.stimContrast(trialData.conditionnumber==conds(k)));
      Deflections(k) = mean(trialData.stimDeflection(trialData.conditionnumber==conds(k)));
end
%Change %s so that they dont equal 0 or 100 and ruin the dprime
HitRight(HitRight==0) = 1./(TotalTrials(HitRight==0)+2);
HitRight(HitRight==1) = 1-1./(TotalTrials(HitRight==1)+2);
HitLeft(HitLeft==0) = 1./(TotalTrials(HitLeft==0)+2);
HitLeft(HitLeft==1) = 1-1./(TotalTrials(HitLeft==1)+2);
%Compute dprimes per condition. 
DPR = [];
lefttrials = [1 5 6 7 10 11 14 15 40 42];
%farincond uses normal absolute condition number as index !
farincond   = [repmat(1,1,15) repmat(16,1,21-16+1) 1 1 1 1 16 16 16 1 1 1 1 1 1 1 1 16 16 16 16 16 16];
for k = lefttrials(ismember(lefttrials,usedconds))
    DPR(k) = norminv(HitLeft(k)) - norminv(HitLeft(farincond(k))) ;
end
%Compute reduction in dprime per pairs of laser ON/OFF
DRED = [];
if strcmp(formatt,'Rec2LR')
DRED(1) = ( DPR(5)-DPR(40) )/ DPR(5); %LOW+FULL
else
DRED(1) = nan; %LOW+FULL 
end
DRED(2) = nan; %Never tested

%Compute average RT
RELRT_DRED = [];
usedconds = unique(trialData.conditionnumber);
RT = cell(max(usedconds),2);
for k=1:length(trialData.stimStart)
    iCond = trialData.conditionnumber(k);
    LTIW = trialData.lickTime{k}(trialData.lickTime{k}>=trialData.respwinStart(k) & trialData.lickTime{k}<=(trialData.respwinEnd(k)+0.002));
    if ~isempty(LTIW) && iCond~=0 
        if trialData.correctResponse(k)
        RT{iCond,1} = [ RT{iCond,1} (LTIW(1) - trialData.stimStart(k)) ];
        else %elseif trialData.firstIncorrect(k)
        RT{iCond,2} = [ RT{iCond,2} (LTIW(1) - trialData.stimStart(k)) ];    
        end
    end
end
clear ans iCond k LTIW usedconds
RELRT_DRED(1) = 0.240       - mean(RT{5,1}); %delayed-low
RELRT_DRED(2) = 0.240       - mean(RT{6,1}); %delayed-high
%RELRT_DRED(3) = 0      - mean(RT{12,1}); %delayed-low MULTI
%RELRT_DRED(4) = 0      - mean(RT{13,1}); %delayed-high MULTI

%Fill BehLatencies struct per condition
BehLat(cc).Task             = 'tMulti';       %tUni/tMulti                
BehLat(cc).Saliency         = 'Low';        %Low/High
BehLat(cc).Silencing        = 'Late';      %Off/Early/Late
BehLat(cc).Modality         = 'Visual';
BehLat(cc).Side             = 'Ipsi'; 
BehLat(cc).LaserOnset       = 0.240;            %Nan/0/0.240 (in s)
BehLat(cc).DprimeCtrl       = DPR(5);
if strcmp(formatt,'Rec2LR')
BehLat(cc).Dprime           = DPR(40);
BehLat(cc).DprimeReduction  = ( DPR(5)-DPR(40) )/ DPR(5);
BehLat(cc).RTmean           = mean(RT{40,1});
else
BehLat(cc).Dprime           = nan;
BehLat(cc).DprimeReduction  = nan;
BehLat(cc).RTmean           = nan;
end
BehLat(cc).RTctrl           = mean(RT{5,1});
BehLat(cc).RTctrlmed        = median(RT{5,1});
BehLat(cc).SessionName      = files_delayed(iSess).name;
%
BehLat(cc+1).Task             = 'tMulti';       %tUni/tMulti                
BehLat(cc+1).Saliency         = 'High';        %Low/High
BehLat(cc+1).Silencing        = 'Late';      %Off/Early/Late
BehLat(cc+1).Modality         = 'Visual';
BehLat(cc+1).Side             = 'Ipsi'; 
BehLat(cc+1).LaserOnset       = 0.240;            %Nan/0/0.240 (in s)
BehLat(cc+1).Dprime           = nan;
BehLat(cc+1).DprimeCtrl       = DPR(6);
BehLat(cc+1).DprimeReduction  = nan;
BehLat(cc+1).RTmean           = nan;
BehLat(cc+1).RTctrl           = mean(RT{6,1});
BehLat(cc+1).RTctrlmed        = median(RT{6,1});
BehLat(cc+1).SessionName      = files_delayed(iSess).name;
%
cc= cc+2;
end
clearvars -except BehLat

%% Get dprimes for MST Tactile
folder_full = 'C:\Scratch\Beh\Data\V1 second bump silencing\V1 full stim silencing';
folder_delayed = 'C:\Scratch\Beh\Data\V1 second bump silencing\V1 delayed silencing';
files_full = dir(strcat(folder_full,'\*.mat'));
files_delayed = dir(strcat(folder_delayed,'\*.mat'));
cc= length(BehLat)+1 ;
%For full silencing
for iSess=1:length(files_full)
clearvars -except iSess files* fi* ALL* cc BehLat
%Load session trialData (no need for neuronal activity)
load(fullfile(files_full(iSess).folder,files_full(iSess).name))

%Reparse and pre-process
pie_beh_reparse
trialData.stimDeflection(isnan(trialData.stimDeflection))=0;
trialData.stimContrast(isnan(trialData.stimContrast))=0;
if isfield(trialData,'stimVISSide')
    trialData.stimDeflection(cell2mat(trialData.stimTACSide)=='L') =  trialData.stimDeflection(cell2mat(trialData.stimTACSide)=='L')*-1;
    trialData.stimContrast(cell2mat(trialData.stimVISSide)=='L')   =  trialData.stimContrast(cell2mat(trialData.stimVISSide)=='L')*-1;    
else
    trialData.stimDeflection(trialData.leftCorrect==1) =  trialData.stimDeflection(trialData.leftCorrect==1)*-1;
    trialData.stimContrast(trialData.leftCorrect==1)   =  trialData.stimContrast(trialData.leftCorrect==1)*-1;
end
%Compute variables per condition
conds = 1:42; %unique(trialData.conditionnumber);
condsnummax = max(conds);
usedconds = unique(trialData.conditionnumber);
for k = 1:condsnummax
      TotalTrials(k) = nansum(trialData.conditionnumber==conds(k));  
      HitRight(k) = nansum(trialData.firstRespLFR=='R' & trialData.conditionnumber==conds(k)) / TotalTrials(k);
      HitLeft(k) = nansum(trialData.firstRespLFR=='L' & trialData.conditionnumber==conds(k)) / TotalTrials(k);
      NoResp(k)  = nansum(isnan(trialData.firstRespLFR) & trialData.conditionnumber==conds(k)) / TotalTrials(k);
      Contrasts(k) = mean(trialData.stimContrast(trialData.conditionnumber==conds(k)));
      Deflections(k) = mean(trialData.stimDeflection(trialData.conditionnumber==conds(k)));
      %CIM
      CorrectCh(k) = nansum(trialData.correctResponse==1 & trialData.conditionnumber==conds(k)) / TotalTrials(k);
      IncorrectCh(k) = nansum(trialData.firstIncorrect==1 & trialData.conditionnumber==conds(k)) / TotalTrials(k);
      NoRespCh(k) = nansum(trialData.noResponse==1 & trialData.conditionnumber==conds(k)) / TotalTrials(k);
end
%Change %s so that they dont equal 0 or 100 and ruin the dprime
HitRight(HitRight==0) = 1./(TotalTrials(HitRight==0)+2);
HitRight(HitRight==1) = 1-1./(TotalTrials(HitRight==1)+2);
HitLeft(HitLeft==0) = 1./(TotalTrials(HitLeft==0)+2);
HitLeft(HitLeft==1) = 1-1./(TotalTrials(HitLeft==1)+2);
%Compute dprimes per condition. 
DPR = nan(1,42);
righttrials = [1  2  3  4  8  9 12 13 16 17 18 19 20 21 29 30 31 32 37 38 39 40 41 42];
%farincond uses normal absolute condition number as index !
farincond   = [repmat(1,1,15) repmat(16,1,21-16+1) 1 1 1 1 16 16 16 1 1 1 1 1 1 1 1 16 16 16 16 16 16];
for k = righttrials(ismember(righttrials,usedconds))
    DPR(k) = norminv(HitRight(k)) - norminv(HitRight(farincond(k))) ;
end
%Compute reduction in dprime per pairs of laser ON/OFF
DRED = [];
DRED(1) = ( DPR(8)-DPR(19) )/ DPR(8); %LOW+FULL
DRED(2) = ( DPR(9)-DPR(38) )/ DPR(9); %HIGH+FULL
% DRED(3) = ( DPR(2)-DPR(27) )/ DPR(2); %LOW+FULL MULTI
% DRED(4) = ( DPR(3)-DPR(28) )/ DPR(3); %HIGH+FULL  MULTI

%Compute average RT
RELRT_DRED = [];
usedconds = unique(trialData.conditionnumber);
RT = cell(42,2);
for k=1:length(trialData.stimStart)
    iCond = trialData.conditionnumber(k);
    LTIW = trialData.lickTime{k}(trialData.lickTime{k}>=trialData.respwinStart(k) & trialData.lickTime{k}<=(trialData.respwinEnd(k)+0.002));
    if ~isempty(LTIW) && iCond~=0 
        if trialData.correctResponse(k)
        RT{iCond,1} = [ RT{iCond,1} (LTIW(1) - trialData.stimStart(k)) ];
        else %elseif trialData.firstIncorrect(k)
        RT{iCond,2} = [ RT{iCond,2} (LTIW(1) - trialData.stimStart(k)) ];    
        end
    end
end
clear ans iCond k LTIW usedconds
RELRT_DRED(1) = 0       - mean(RT{8,1}); %full-low
RELRT_DRED(2) = 0       - mean(RT{9,1}); %full-high
%RELRT_DRED(3) = 0      - mean(RT{12,1}); %full-low MULTI
%RELRT_DRED(4) = 0      - mean(RT{13,1}); %full-high MULTI

%Fill BehLatencies struct per condition
BehLat(cc).Task             = 'tMulti';       %tUni/tMulti                
BehLat(cc).Saliency         = 'Low';        %Low/High
BehLat(cc).Silencing        = 'Early';      %Off/Early/Late
BehLat(cc).Modality         = 'Tactile';      
BehLat(cc).Side             = 'Contra';     
BehLat(cc).LaserOnset       = 0;            %Nan/0/0.240 (in s)
BehLat(cc).Dprime           = DPR(19);
BehLat(cc).DprimeCtrl       = DPR(8);
BehLat(cc).DprimeReduction  = ( DPR(8)-DPR(19) )/ DPR(8);
BehLat(cc).RTmean           = mean(RT{19,1});
BehLat(cc).RTctrl           = mean(RT{8,1});
BehLat(cc).RTctrlmed        = median(RT{8,1});
BehLat(cc).SessionName      = files_full(iSess).name;
%
BehLat(cc+1).Task             = 'tMulti';       %tUni/tMulti                
BehLat(cc+1).Saliency         = 'High';        %Low/High
BehLat(cc+1).Silencing        = 'Early';      %Off/Early/Late
BehLat(cc+1).Modality         = 'Tactile';      
BehLat(cc+1).Side             = 'Contra';  
BehLat(cc+1).LaserOnset       = 0;            %Nan/0/0.240 (in s)
BehLat(cc+1).Dprime           = DPR(38);
BehLat(cc+1).DprimeCtrl       = DPR(9);
BehLat(cc+1).DprimeReduction  = ( DPR(9)-DPR(38) )/ DPR(9);
BehLat(cc+1).RTmean           = mean(RT{38,1});
BehLat(cc+1).RTctrl           = mean(RT{9,1});
BehLat(cc+1).RTctrlmed        = median(RT{9,1});
BehLat(cc+1).SessionName      = files_full(iSess).name;
%
BehLat(cc+2).Task             = 'tMulti';       %tUni/tMulti                
BehLat(cc+2).Saliency         = 'Low';        %Low/High
BehLat(cc+2).Silencing        = 'Early';      %Off/Early/Late
BehLat(cc+2).Modality         = 'Tactile';      
BehLat(cc+2).Side             = 'Ipsi';     
BehLat(cc+2).LaserOnset       = 0;            %Nan/0/0.240 (in s)
BehLat(cc+2).Dprime           = DPR(26);
BehLat(cc+2).DprimeCtrl       = DPR(10);
BehLat(cc+2).DprimeReduction  = ( DPR(10)-DPR(26) )/ DPR(10);
BehLat(cc+2).RTmean           = mean(RT{26,1});
BehLat(cc+2).RTctrl           = mean(RT{10,1});
BehLat(cc+2).RTctrlmed        = median(RT{10,1});
BehLat(cc+2).SessionName      = files_full(iSess).name;
%
BehLat(cc+3).Task             = 'tMulti';       %tUni/tMulti                
BehLat(cc+3).Saliency         = 'High';        %Low/High
BehLat(cc+3).Silencing        = 'Early';      %Off/Early/Late
BehLat(cc+3).Modality         = 'Tactile';      
BehLat(cc+3).Side             = 'Ipsi';  
BehLat(cc+3).LaserOnset       = 0;            %Nan/0/0.240 (in s)
BehLat(cc+3).Dprime           = nan; %DPR(51)
BehLat(cc+3).DprimeCtrl       = DPR(11);
BehLat(cc+3).DprimeReduction  = nan; %( DPR(11)-DPR(51) )/ DPR(11);
BehLat(cc+3).RTmean           = nan; %mean(RT{51,1});
BehLat(cc+3).RTctrl           = mean(RT{11,1});
BehLat(cc+3).RTctrlmed        = median(RT{11,1});
BehLat(cc+3).SessionName      = files_full(iSess).name;
%
cc= cc+4;
end
%For delayed silencing
for iSess=1:length(files_delayed)
clearvars -except iSess files* fi* ALL* cc BehLat
%Load session trialData (no need for neuronal activity)
load(fullfile(files_delayed(iSess).folder,files_delayed(iSess).name))

%Reparse and pre-process
pie_beh_reparse
trialData.stimDeflection(isnan(trialData.stimDeflection))=0;
trialData.stimContrast(isnan(trialData.stimContrast))=0;
if isfield(trialData,'stimVISSide')
    trialData.stimDeflection(cell2mat(trialData.stimTACSide)=='L') =  trialData.stimDeflection(cell2mat(trialData.stimTACSide)=='L')*-1;
    trialData.stimContrast(cell2mat(trialData.stimVISSide)=='L')   =  trialData.stimContrast(cell2mat(trialData.stimVISSide)=='L')*-1;    
else
    trialData.stimDeflection(trialData.leftCorrect==1) =  trialData.stimDeflection(trialData.leftCorrect==1)*-1;
    trialData.stimContrast(trialData.leftCorrect==1)   =  trialData.stimContrast(trialData.leftCorrect==1)*-1;
end
%Compute variables per condition
conds = 1:42; %unique(trialData.conditionnumber);
condsnummax = max(conds);
usedconds = unique(trialData.conditionnumber);
for k = 1:condsnummax
      TotalTrials(k) = nansum(trialData.conditionnumber==conds(k));  
      HitRight(k) = nansum(trialData.firstRespLFR=='R' & trialData.conditionnumber==conds(k)) / TotalTrials(k);
      HitLeft(k) = nansum(trialData.firstRespLFR=='L' & trialData.conditionnumber==conds(k)) / TotalTrials(k);
      NoResp(k)  = nansum(isnan(trialData.firstRespLFR) & trialData.conditionnumber==conds(k)) / TotalTrials(k);
      Contrasts(k) = mean(trialData.stimContrast(trialData.conditionnumber==conds(k)));
      Deflections(k) = mean(trialData.stimDeflection(trialData.conditionnumber==conds(k)));
      %CIM
      CorrectCh(k) = nansum(trialData.correctResponse==1 & trialData.conditionnumber==conds(k)) / TotalTrials(k);
      IncorrectCh(k) = nansum(trialData.firstIncorrect==1 & trialData.conditionnumber==conds(k)) / TotalTrials(k);
      NoRespCh(k) = nansum(trialData.noResponse==1 & trialData.conditionnumber==conds(k)) / TotalTrials(k);
end
%Change %s so that they dont equal 0 or 100 and ruin the dprime
HitRight(HitRight==0) = 1./(TotalTrials(HitRight==0)+2);
HitRight(HitRight==1) = 1-1./(TotalTrials(HitRight==1)+2);
HitLeft(HitLeft==0) = 1./(TotalTrials(HitLeft==0)+2);
HitLeft(HitLeft==1) = 1-1./(TotalTrials(HitLeft==1)+2);
%Compute dprimes per condition. 
DPR = nan(1,42);
righttrials = [1  2  3  4  8  9 12 13 16 17 18 19 20 21 29 30 31 32 37 38 39 40 41 42];
%farincond uses normal absolute condition number as index !
farincond   = [repmat(1,1,15) repmat(16,1,21-16+1) 1 1 1 1 16 16 16 1 1 1 1 1 1 1 1 16 16 16 16 16 16];
for k = righttrials(ismember(righttrials,usedconds))
    DPR(k) = norminv(HitRight(k)) - norminv(HitRight(farincond(k))) ;
end
%Compute reduction in dprime per pairs of laser ON/OFF
DRED = [];
DRED(1) = ( DPR(8)-DPR(19) )/ DPR(8); %LOW+FULL
DRED(2) = ( DPR(9)-DPR(38) )/ DPR(9); %HIGH+FULL
% DRED(3) = ( DPR(2)-DPR(27) )/ DPR(2); %LOW+FULL MULTI
% DRED(4) = ( DPR(3)-DPR(28) )/ DPR(3); %HIGH+FULL  MULTI

%Compute average RT
RELRT_DRED = [];
usedconds = unique(trialData.conditionnumber);
RT = cell(42,2);
for k=1:length(trialData.stimStart)
    iCond = trialData.conditionnumber(k);
    LTIW = trialData.lickTime{k}(trialData.lickTime{k}>=trialData.respwinStart(k) & trialData.lickTime{k}<=(trialData.respwinEnd(k)+0.002));
    if ~isempty(LTIW) && iCond~=0 
        if trialData.correctResponse(k)
        RT{iCond,1} = [ RT{iCond,1} (LTIW(1) - trialData.stimStart(k)) ];
        else %elseif trialData.firstIncorrect(k)
        RT{iCond,2} = [ RT{iCond,2} (LTIW(1) - trialData.stimStart(k)) ];    
        end
    end
end
clear ans iCond k LTIW usedconds
RELRT_DRED(1) = 0.240       - mean(RT{8,1}); %delayed-low
RELRT_DRED(2) = 0.240       - mean(RT{9,1}); %delayed-high
%RELRT_DRED(3) = 0      - mean(RT{12,1}); %delayed-low MULTI
%RELRT_DRED(4) = 0      - mean(RT{13,1}); %delayed-high MULTI

%Fill BehLatencies struct per condition
BehLat(cc).Task             = 'tMulti';       %tUni/tMulti                
BehLat(cc).Saliency         = 'Low';        %Low/High
BehLat(cc).Silencing        = 'Late';      %Off/Early/Late
BehLat(cc).Modality         = 'Tactile';      
BehLat(cc).Side             = 'Contra';  
BehLat(cc).LaserOnset       = 0.240;            %Nan/0/0.240 (in s)
BehLat(cc).Dprime           = DPR(19);
BehLat(cc).DprimeCtrl       = DPR(8);
BehLat(cc).DprimeReduction  = ( DPR(8)-DPR(19) )/ DPR(8);
BehLat(cc).RTmean           = mean(RT{19,1});
BehLat(cc).RTctrl           = mean(RT{8,1});
BehLat(cc).RTctrlmed        = median(RT{8,1});
BehLat(cc).SessionName      = files_delayed(iSess).name;
%
BehLat(cc+1).Task             = 'tMulti';       %tUni/tMulti                
BehLat(cc+1).Saliency         = 'High';        %Low/High
BehLat(cc+1).Silencing        = 'Late';      %Off/Early/Late
BehLat(cc+1).Modality         = 'Tactile';      
BehLat(cc+1).Side             = 'Contra';  
BehLat(cc+1).LaserOnset       = 0.240;            %Nan/0/0.240 (in s)
BehLat(cc+1).Dprime           = DPR(38);
BehLat(cc+1).DprimeCtrl       = DPR(9);
BehLat(cc+1).DprimeReduction  = ( DPR(9)-DPR(38) )/ DPR(9);
BehLat(cc+1).RTmean           = mean(RT{38,1});
BehLat(cc+1).RTctrl           = mean(RT{9,1});
BehLat(cc+1).RTctrlmed        = median(RT{9,1});
BehLat(cc+1).SessionName      = files_delayed(iSess).name;
%
BehLat(cc+2).Task             = 'tMulti';       %tUni/tMulti                
BehLat(cc+2).Saliency         = 'Low';        %Low/High
BehLat(cc+2).Silencing        = 'Late';      %Off/Early/Late
BehLat(cc+2).Modality         = 'Tactile';      
BehLat(cc+2).Side             = 'Ipsi';     
BehLat(cc+2).LaserOnset       = 0.240;            %Nan/0/0.240 (in s)
BehLat(cc+2).Dprime           = DPR(26);
BehLat(cc+2).DprimeCtrl       = DPR(10);
BehLat(cc+2).DprimeReduction  = ( DPR(10)-DPR(26) )/ DPR(10);
BehLat(cc+2).RTmean           = mean(RT{26,1});
BehLat(cc+2).RTctrl           = mean(RT{10,1});
BehLat(cc+2).RTctrlmed        = median(RT{10,1});
BehLat(cc+2).SessionName      = files_delayed(iSess).name;
%
BehLat(cc+3).Task             = 'tMulti';       %tUni/tMulti                
BehLat(cc+3).Saliency         = 'High';        %Low/High
BehLat(cc+3).Silencing        = 'Late';      %Off/Early/Late
BehLat(cc+3).Modality         = 'Tactile';      
BehLat(cc+3).Side             = 'Ipsi';  
BehLat(cc+3).LaserOnset       = 0.240;            %Nan/0/0.240 (in s)
BehLat(cc+3).Dprime           = nan; %DPR(51);
BehLat(cc+3).DprimeCtrl       = DPR(11);
BehLat(cc+3).DprimeReduction  = nan; %( DPR(11)-DPR(51) )/ DPR(11);
BehLat(cc+3).RTmean           = nan; %mean(RT{51,1});
BehLat(cc+3).RTctrl           = mean(RT{11,1});
BehLat(cc+3).RTctrlmed        = median(RT{11,1});
BehLat(cc+3).SessionName      = files_delayed(iSess).name;
%
cc= cc+4;
end
clearvars -except BehLat

%% Replace the 0 dprimes that should be nans (no max cond in)
BehLat([30]).Dprime = nan;
BehLat([32]).Dprime = nan;
BehLat([42]).Dprime = nan;
BehLat([44]).Dprime = nan;