%This code generates Supplementary Figure S8c from oude Lohuis et al., 2022.
%Reaction times are obtained on a trial-by-trial basis for different
%conditions. Plotted as bar-plots. LME stats.

%% Insert directory where OudeLohuisetal_2022_NatComms is saved:
maindirectory = 'C:\Users\jeanl\surfdrive\Shared\Manuscript - 2nd Bump\Submission - NatComm\rev2';

%% Load BehLat

load( strcat(maindirectory,'\OudeLohuisetal_2022_NatComms\6VisuotactileTask\Data6_3\BehLat.mat') );
%prep ------------------------------------------------------------------
sessionnames = categorical({BehLat.SessionName});
sessionnum = double(sessionnames);
for i=1:length(BehLat)
    SN = BehLat(i).SessionName;
    BehLat(i).Mouse = SN(strfind(SN,'J'):strfind(SN,'J')+3);
end

%% Plot
%For UST only use early since early and late same session
UV_Thr = [BehLat(strcmp({BehLat.Task},'tUni') & strcmp({BehLat.Saliency},'Low') & strcmp({BehLat.Modality},'Visual') & strcmp({BehLat.Silencing},'Early')).RTctrlmed];
UV_Max = [BehLat(strcmp({BehLat.Task},'tUni') & strcmp({BehLat.Saliency},'High') & strcmp({BehLat.Modality},'Visual') & strcmp({BehLat.Silencing},'Early')).RTctrlmed];
MV_Thr = [BehLat(strcmp({BehLat.Task},'tMulti') & strcmp({BehLat.Saliency},'Low') & strcmp({BehLat.Modality},'Visual') ).RTctrlmed];
MV_Max = [BehLat(strcmp({BehLat.Task},'tMulti') & strcmp({BehLat.Saliency},'High') & strcmp({BehLat.Modality},'Visual') ).RTctrlmed];
MT_Thr = [BehLat(strcmp({BehLat.Task},'tMulti') & strcmp({BehLat.Saliency},'Low') & strcmp({BehLat.Modality},'Tactile') ).RTctrlmed];
MT_Max = [BehLat(strcmp({BehLat.Task},'tMulti') & strcmp({BehLat.Saliency},'High') & strcmp({BehLat.Modality},'Tactile') ).RTctrlmed];

% figure(); hold on
% bar(1,nanmedian(UV_Thr))
% bar(2,nanmedian(UV_Max))
% bar(3,nanmedian(MV_Thr))
% bar(4,nanmedian(MV_Max))
% bar(5,nanmedian(MT_Thr))
% bar(6,nanmedian(MT_Max))

%% Boxplots
data = [[UV_Thr';nan(16,1)] [UV_Max';nan(16,1)] MV_Thr' MV_Max' MT_Thr' MT_Max'];
figure(); hold on
boxplot(data,'plotstyle','compact','Color',[0 0 .4; 0 0 .9; 0 0 .4; 0 0 .9; 0 .4 0; 0 .9 0],'boxstyle','filled', 'medianstyle','target','whisker',1.5,'outliersize',0.0001,'Widths',5)
ylim([0.2 1])

%% LME Stats
%LME stats are done using single trials, as follows:

%% Get Reaction Times for: Visual Contralateral
%% UST
folder = strcat(maindirectory,'\OudeLohuisetal_2022_NatComms\6VisuotactileTask\Data6_3\UST');
files = dir(strcat(folder,'\*.mat'));
tbl = table;

for iSess=1:length(files)
clearvars -except iSess files fi* fo* tbl* maindirectory
%Load session trialData (no need for neuronal activity)
load(fullfile(files(iSess).folder,files(iSess).name))
pie_beh_reparse
%Get mouse name
nnn = strfind(files(iSess).name,'_J');
mousename = files(iSess).name(nnn+1:nnn+4);
%Get RTs
usedconds = unique(trialData.conditionnumber);
sessRTs = cell(max(usedconds),2);
for k=1:length(trialData.stimStart)
    iCond = trialData.conditionnumber(k);
    LTIW = trialData.lickTime{k}(trialData.lickTime{k}>=trialData.stimStart(k));
    if ~isempty(LTIW) && iCond~=0 
         if (LTIW(1) - trialData.stimStart(k))>=0.2 && (LTIW(1) - trialData.stimStart(k))<=1 %Remove licks before 200ms - (too fast for a proper task-engaged trial)
            if trialData.correctResponse(k)
            sessRTs{iCond,1} = [ sessRTs{iCond,1} (LTIW(1) - trialData.stimStart(k)) ];
            elseif trialData.firstIncorrect(k)
            sessRTs{iCond,2} = [ sessRTs{iCond,2} (LTIW(1) - trialData.stimStart(k)) ];    
            end
         end
    end
end
%Conditions needed:
%VIS Contra Thr : 2
tbl2 = table;
tbl2.RT = sessRTs{2,1}';
tbl2.Mouse      = repmat({mousename},size(tbl2,1),1);
tbl2.Session    = repmat(iSess,size(tbl2,1),1);
tbl2.Saliency   = repmat({'Thr'},size(tbl2,1),1);
tbl2.Side       = repmat({'Contra'},size(tbl2,1),1);
tbl2.Task       = repmat({'UST'},size(tbl2,1),1);
tbl = [tbl;tbl2];
%VIS Contra Max : 3
tbl2 = table;
tbl2.RT = sessRTs{3,1}';
tbl2.Mouse      = repmat({mousename},size(tbl2,1),1);
tbl2.Session    = repmat(iSess,size(tbl2,1),1);
tbl2.Saliency   = repmat({'Max'},size(tbl2,1),1);
tbl2.Side       = repmat({'Contra'},size(tbl2,1),1);
tbl2.Task       = repmat({'UST'},size(tbl2,1),1);
tbl = [tbl;tbl2];
%VIS Ipsi Thr   : 5
tbl2 = table;
tbl2.RT = sessRTs{5,1}';
tbl2.Mouse      = repmat({mousename},size(tbl2,1),1);
tbl2.Session    = repmat(iSess,size(tbl2,1),1);
tbl2.Saliency   = repmat({'Thr'},size(tbl2,1),1);
tbl2.Side       = repmat({'Ipsi'},size(tbl2,1),1);
tbl2.Task       = repmat({'UST'},size(tbl2,1),1);
tbl = [tbl;tbl2];
%VIS Ipsi Max   : 6
tbl2 = table;
tbl2.RT = sessRTs{6,1}';
tbl2.Mouse      = repmat({mousename},size(tbl2,1),1);
tbl2.Session    = repmat(iSess,size(tbl2,1),1);
tbl2.Saliency   = repmat({'Max'},size(tbl2,1),1);
tbl2.Side       = repmat({'Ipsi'},size(tbl2,1),1);
tbl2.Task       = repmat({'UST'},size(tbl2,1),1);
tbl = [tbl;tbl2];
%TAC isn't present for UST
end
tblUST = tbl;
clearvars tbl tbl2 nnn sessRTs 

%% %%%%%%%Contralateral MST
folder_full = strcat(maindirectory,'\OudeLohuisetal_2022_NatComms\6VisuotactileTask\Data6_3\MST\V1 full stim silencing');
folder_delayed = strcat(maindirectory,'\OudeLohuisetal_2022_NatComms\6VisuotactileTask\Data6_3\MST\V1 delayed silencing');
files_full = dir(strcat(folder_full,'\*.mat'));
files_delayed = dir(strcat(folder_delayed,'\*.mat'));

files = files_full;
folder = folder_full;
tbl = table;
for iSess=1:length(files)
clearvars -except iSess files fi* fo* tbl* maindirectory
%Load session trialData (no need for neuronal activity)
load(fullfile(files(iSess).folder,files(iSess).name))
pie_beh_reparse
%Get mouse name
nnn = strfind(files(iSess).name,'_J');
mousename = files(iSess).name(nnn+1:nnn+4);
%Get RTs
usedconds = unique(trialData.conditionnumber);
sessRTs = cell(max(usedconds),2);
for k=1:length(trialData.stimStart)
    iCond = trialData.conditionnumber(k);
    LTIW = trialData.lickTime{k}(trialData.lickTime{k}>=trialData.stimStart(k));
    if ~isempty(LTIW) && iCond~=0 
         if (LTIW(1) - trialData.stimStart(k))>=0.2 && (LTIW(1) - trialData.stimStart(k))<=1 %Remove licks before 200ms - (too fast for a proper task-engaged trial)
            if trialData.correctResponse(k)
            sessRTs{iCond,1} = [ sessRTs{iCond,1} (LTIW(1) - trialData.stimStart(k)) ];
            elseif trialData.firstIncorrect(k)
            sessRTs{iCond,2} = [ sessRTs{iCond,2} (LTIW(1) - trialData.stimStart(k)) ];    
            end
         end
    end
end
%Conditions needed:
%VIS Contra Thr : 2
tbl2 = table;
tbl2.RT = sessRTs{2,1}';
tbl2.Mouse      = repmat({mousename},size(tbl2,1),1);
tbl2.Session    = repmat(iSess,size(tbl2,1),1);
tbl2.Saliency   = repmat({'Thr'},size(tbl2,1),1);
tbl2.Side       = repmat({'Contra'},size(tbl2,1),1);
tbl2.Task       = repmat({'MST'},size(tbl2,1),1);
tbl = [tbl;tbl2];
%VIS Contra Max : 3
tbl2 = table;
tbl2.RT = sessRTs{3,1}';
tbl2.Mouse      = repmat({mousename},size(tbl2,1),1);
tbl2.Session    = repmat(iSess,size(tbl2,1),1);
tbl2.Saliency   = repmat({'Max'},size(tbl2,1),1);
tbl2.Side       = repmat({'Contra'},size(tbl2,1),1);
tbl2.Task       = repmat({'MST'},size(tbl2,1),1);
tbl = [tbl;tbl2];
%VIS Ipsi Thr   : 5
tbl2 = table;
tbl2.RT = sessRTs{5,1}';
tbl2.Mouse      = repmat({mousename},size(tbl2,1),1);
tbl2.Session    = repmat(iSess,size(tbl2,1),1);
tbl2.Saliency   = repmat({'Thr'},size(tbl2,1),1);
tbl2.Side       = repmat({'Ipsi'},size(tbl2,1),1);
tbl2.Task       = repmat({'MST'},size(tbl2,1),1);
tbl = [tbl;tbl2];
%VIS Ipsi Max   : 6
tbl2 = table;
tbl2.RT = sessRTs{6,1}';
tbl2.Mouse      = repmat({mousename},size(tbl2,1),1);
tbl2.Session    = repmat(iSess,size(tbl2,1),1);
tbl2.Saliency   = repmat({'Max'},size(tbl2,1),1);
tbl2.Side       = repmat({'Ipsi'},size(tbl2,1),1);
tbl2.Task       = repmat({'MST'},size(tbl2,1),1);
tbl = [tbl;tbl2];
%TAC isn't present for UST
end
tblMST1 = tbl;
clearvars tbl tbl2 nnn sessRTs

files = files_delayed;
folder = folder_delayed;
tbl = table;
for iSess=1:length(files)
clearvars -except iSess files fi* fo* tbl* maindirectory
%Load session trialData (no need for neuronal activity)
load(fullfile(files(iSess).folder,files(iSess).name))
pie_beh_reparse
%Get mouse name
nnn = strfind(files(iSess).name,'_J');
mousename = files(iSess).name(nnn+1:nnn+4);
%Get RTs
usedconds = unique(trialData.conditionnumber);
sessRTs = cell(max(usedconds),2);
for k=1:length(trialData.stimStart)
    iCond = trialData.conditionnumber(k);
    LTIW = trialData.lickTime{k}(trialData.lickTime{k}>=trialData.stimStart(k));
    if ~isempty(LTIW) && iCond~=0 
         if (LTIW(1) - trialData.stimStart(k))>=0.2 && (LTIW(1) - trialData.stimStart(k))<=1 %Remove licks before 200ms - (too fast for a proper task-engaged trial)
            if trialData.correctResponse(k)
            sessRTs{iCond,1} = [ sessRTs{iCond,1} (LTIW(1) - trialData.stimStart(k)) ];
            elseif trialData.firstIncorrect(k)
            sessRTs{iCond,2} = [ sessRTs{iCond,2} (LTIW(1) - trialData.stimStart(k)) ];    
            end
         end
    end
end
%Conditions needed:
%VIS Contra Thr : 2
tbl2 = table;
tbl2.RT = sessRTs{2,1}';
tbl2.Mouse      = repmat({mousename},size(tbl2,1),1);
tbl2.Session    = repmat(iSess,size(tbl2,1),1);
tbl2.Saliency   = repmat({'Thr'},size(tbl2,1),1);
tbl2.Side       = repmat({'Contra'},size(tbl2,1),1);
tbl2.Task       = repmat({'MST'},size(tbl2,1),1);
tbl = [tbl;tbl2];
%VIS Contra Max : 3
tbl2 = table;
tbl2.RT = sessRTs{3,1}';
tbl2.Mouse      = repmat({mousename},size(tbl2,1),1);
tbl2.Session    = repmat(iSess,size(tbl2,1),1);
tbl2.Saliency   = repmat({'Max'},size(tbl2,1),1);
tbl2.Side       = repmat({'Contra'},size(tbl2,1),1);
tbl2.Task       = repmat({'MST'},size(tbl2,1),1);
tbl = [tbl;tbl2];
%VIS Ipsi Thr   : 5
tbl2 = table;
tbl2.RT = sessRTs{5,1}';
tbl2.Mouse      = repmat({mousename},size(tbl2,1),1);
tbl2.Session    = repmat(iSess,size(tbl2,1),1);
tbl2.Saliency   = repmat({'Thr'},size(tbl2,1),1);
tbl2.Side       = repmat({'Ipsi'},size(tbl2,1),1);
tbl2.Task       = repmat({'MST'},size(tbl2,1),1);
tbl = [tbl;tbl2];
%VIS Ipsi Max   : 6
tbl2 = table;
tbl2.RT = sessRTs{6,1}';
tbl2.Mouse      = repmat({mousename},size(tbl2,1),1);
tbl2.Session    = repmat(iSess,size(tbl2,1),1);
tbl2.Saliency   = repmat({'Max'},size(tbl2,1),1);
tbl2.Side       = repmat({'Ipsi'},size(tbl2,1),1);
tbl2.Task       = repmat({'MST'},size(tbl2,1),1);
tbl = [tbl;tbl2];
%TAC isn't present for UST
end
tblMST2 = tbl;
clearvars tbl tbl2 nnn sessRTs 

%% Pool together
tbl = [tblUST; tblMST1; tblMST2];
tbl.Mouse = categorical(tbl.Mouse);
tbl.Session = categorical(tbl.Session);
tbl.Saliency = categorical(tbl.Saliency);
tbl.Side = categorical(tbl.Side);
tbl.Task = categorical(tbl.Task);


%% STATS
% lme     = fitglme(tbl,'RT ~ Task + Saliency + Side + (1|Session:Mouse)',...
lme     = fitglme(tbl,'RT ~ Task + Saliency + Side',...
                        'Distribution','InverseGaussian','Link','identity');
stats   = dataset2table(anova(lme,'DFMethod','Residual'));
% 
% %% Do stats
% lme = fitlme(tbl,'RT ~ Task + Saliency + Side + (1|Mouse) + (1|Session:Mouse)');
% stats = dataset2table(anova(lme,'DFMethod','Satterthwaite'));
% p1 = stats{2,5};

%% Separate per Cohort to test Thr vs Max
tbl = [tblUST; tblMST1; tblMST2];
D = table2struct(tbl);
D = D(strcmp({D.Task},'UST'));
tbl2 = struct2table(D);
tbl2.Mouse = categorical(tbl2.Mouse);
tbl2.Session = categorical(tbl2.Session);
tbl2.Saliency = categorical(tbl2.Saliency);
tbl2.Side = categorical(tbl2.Side);
tbl2.Task = categorical(tbl2.Task);
lme     = fitglme(tbl2,'RT ~ Saliency + Side + (1|Session:Mouse)',...
                        'Distribution','InverseGaussian','Link','identity');
stats   = dataset2table(anova(lme,'DFMethod','Residual'));

D = table2struct(tbl);
D = D(strcmp({D.Task},'MST'));
tbl2 = struct2table(D);
tbl2.Mouse = categorical(tbl2.Mouse);
tbl2.Session = categorical(tbl2.Session);
tbl2.Saliency = categorical(tbl2.Saliency);
tbl2.Side = categorical(tbl2.Side);
tbl2.Task = categorical(tbl2.Task);
lme     = fitglme(tbl2,'RT ~ Saliency + Side + (1|Session:Mouse)',...
                        'Distribution','InverseGaussian','Link','identity');
stats   = dataset2table(anova(lme,'DFMethod','Residual'));

%% Now tactile
%% %%%%%%% MST
folder_full = strcat(maindirectory,'\OudeLohuisetal_2022_NatComms\6VisuotactileTask\Data6_3\MST\V1 full stim silencing');
folder_delayed = strcat(maindirectory,'\OudeLohuisetal_2022_NatComms\6VisuotactileTask\Data6_3\MST\V1 delayed silencing');
files_full = dir(strcat(folder_full,'\*.mat'));
files_delayed = dir(strcat(folder_delayed,'\*.mat'));

files = files_full;
folder = folder_full;
tbl = table;
for iSess=1:length(files)
clearvars -except iSess files fi* fo* tbl* maindirectory
%Load session trialData (no need for neuronal activity)
load(fullfile(files(iSess).folder,files(iSess).name))
pie_beh_reparse
%Get mouse name
nnn = strfind(files(iSess).name,'_J');
mousename = files(iSess).name(nnn+1:nnn+4);
%Get RTs
usedconds = unique(trialData.conditionnumber);
sessRTs = cell(max(usedconds),2);
for k=1:length(trialData.stimStart)
    iCond = trialData.conditionnumber(k);
    LTIW = trialData.lickTime{k}(trialData.lickTime{k}>=trialData.stimStart(k));
    if ~isempty(LTIW) && iCond~=0 
         if (LTIW(1) - trialData.stimStart(k))>=0.2 && (LTIW(1) - trialData.stimStart(k))<=1 %Remove licks before 200ms - (too fast for a proper task-engaged trial)
            if trialData.correctResponse(k)
            sessRTs{iCond,1} = [ sessRTs{iCond,1} (LTIW(1) - trialData.stimStart(k)) ];
            elseif trialData.firstIncorrect(k)
            sessRTs{iCond,2} = [ sessRTs{iCond,2} (LTIW(1) - trialData.stimStart(k)) ];    
            end
         end
    end
end
%Conditions needed:
%TAC Contra Thr : 8
tbl2 = table;
tbl2.RT = sessRTs{8,1}';
tbl2.Mouse      = repmat({mousename},size(tbl2,1),1);
tbl2.Session    = repmat(iSess,size(tbl2,1),1);
tbl2.Saliency   = repmat({'Thr'},size(tbl2,1),1);
tbl2.Side       = repmat({'Contra'},size(tbl2,1),1);
tbl2.Task       = repmat({'MST'},size(tbl2,1),1);
tbl = [tbl;tbl2];
%VIS Contra Max : 9
tbl2 = table;
tbl2.RT = sessRTs{9,1}';
tbl2.Mouse      = repmat({mousename},size(tbl2,1),1);
tbl2.Session    = repmat(iSess,size(tbl2,1),1);
tbl2.Saliency   = repmat({'Max'},size(tbl2,1),1);
tbl2.Side       = repmat({'Contra'},size(tbl2,1),1);
tbl2.Task       = repmat({'MST'},size(tbl2,1),1);
tbl = [tbl;tbl2];
%VIS Ipsi Thr   : 10
tbl2 = table;
tbl2.RT = sessRTs{10,1}';
tbl2.Mouse      = repmat({mousename},size(tbl2,1),1);
tbl2.Session    = repmat(iSess,size(tbl2,1),1);
tbl2.Saliency   = repmat({'Thr'},size(tbl2,1),1);
tbl2.Side       = repmat({'Ipsi'},size(tbl2,1),1);
tbl2.Task       = repmat({'MST'},size(tbl2,1),1);
tbl = [tbl;tbl2];
%VIS Ipsi Max   : 11
tbl2 = table;
tbl2.RT = sessRTs{11,1}';
tbl2.Mouse      = repmat({mousename},size(tbl2,1),1);
tbl2.Session    = repmat(iSess,size(tbl2,1),1);
tbl2.Saliency   = repmat({'Max'},size(tbl2,1),1);
tbl2.Side       = repmat({'Ipsi'},size(tbl2,1),1);
tbl2.Task       = repmat({'MST'},size(tbl2,1),1);
tbl = [tbl;tbl2];
%TAC isn't present for UST
end
tblMST1 = tbl;
clearvars tbl tbl2 nnn sessRTs 

files = files_delayed;
folder = folder_delayed;
tbl = table;
for iSess=1:length(files)
clearvars -except iSess files fi* fo* tbl* maindirectory
%Load session trialData (no need for neuronal activity)
load(fullfile(files(iSess).folder,files(iSess).name))
pie_beh_reparse
%Get mouse name
nnn = strfind(files(iSess).name,'_J');
mousename = files(iSess).name(nnn+1:nnn+4);
%Get RTs
usedconds = unique(trialData.conditionnumber);
sessRTs = cell(max(usedconds),2);
for k=1:length(trialData.stimStart)
    iCond = trialData.conditionnumber(k);
    LTIW = trialData.lickTime{k}(trialData.lickTime{k}>=trialData.stimStart(k));
    if ~isempty(LTIW) && iCond~=0 
         if (LTIW(1) - trialData.stimStart(k))>=0.2 && (LTIW(1) - trialData.stimStart(k))<=1 %Remove licks before 200ms - (too fast for a proper task-engaged trial)
            if trialData.correctResponse(k)
            sessRTs{iCond,1} = [ sessRTs{iCond,1} (LTIW(1) - trialData.stimStart(k)) ];
            elseif trialData.firstIncorrect(k)
            sessRTs{iCond,2} = [ sessRTs{iCond,2} (LTIW(1) - trialData.stimStart(k)) ];    
            end
         end
    end
end
%Conditions needed:
%VIS Contra Thr : 8
tbl2 = table;
tbl2.RT = sessRTs{8,1}';
tbl2.Mouse      = repmat({mousename},size(tbl2,1),1);
tbl2.Session    = repmat(iSess,size(tbl2,1),1);
tbl2.Saliency   = repmat({'Thr'},size(tbl2,1),1);
tbl2.Side       = repmat({'Contra'},size(tbl2,1),1);
tbl2.Task       = repmat({'MST'},size(tbl2,1),1);
tbl = [tbl;tbl2];
%VIS Contra Max : 9
tbl2 = table;
tbl2.RT = sessRTs{9,1}';
tbl2.Mouse      = repmat({mousename},size(tbl2,1),1);
tbl2.Session    = repmat(iSess,size(tbl2,1),1);
tbl2.Saliency   = repmat({'Max'},size(tbl2,1),1);
tbl2.Side       = repmat({'Contra'},size(tbl2,1),1);
tbl2.Task       = repmat({'MST'},size(tbl2,1),1);
tbl = [tbl;tbl2];
%VIS Ipsi Thr   : 10
tbl2 = table;
tbl2.RT = sessRTs{10,1}';
tbl2.Mouse      = repmat({mousename},size(tbl2,1),1);
tbl2.Session    = repmat(iSess,size(tbl2,1),1);
tbl2.Saliency   = repmat({'Thr'},size(tbl2,1),1);
tbl2.Side       = repmat({'Ipsi'},size(tbl2,1),1);
tbl2.Task       = repmat({'MST'},size(tbl2,1),1);
tbl = [tbl;tbl2];
%VIS Ipsi Max   : 11
tbl2 = table;
tbl2.RT = sessRTs{11,1}';
tbl2.Mouse      = repmat({mousename},size(tbl2,1),1);
tbl2.Session    = repmat(iSess,size(tbl2,1),1);
tbl2.Saliency   = repmat({'Max'},size(tbl2,1),1);
tbl2.Side       = repmat({'Ipsi'},size(tbl2,1),1);
tbl2.Task       = repmat({'MST'},size(tbl2,1),1);
tbl = [tbl;tbl2];
%TAC isn't present for UST
end
tblMST2 = tbl;
clearvars tbl tbl2 nnn sessRTs 
tbl = [tblMST1 ; tblMST2];
tbl.Mouse = categorical(tbl.Mouse);
tbl.Session = categorical(tbl.Session);
tbl.Saliency = categorical(tbl.Saliency);
tbl.Side = categorical(tbl.Side);
tbl.Task = categorical(tbl.Task);
lme     = fitglme(tbl,'RT ~ Saliency + Side + (1|Session:Mouse)',...
                        'Distribution','InverseGaussian','Link','identity');
stats   = dataset2table(anova(lme,'DFMethod','Residual'));

