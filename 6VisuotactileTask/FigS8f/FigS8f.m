%FigS8f 
%Supplementary Figure S8f can be generated using the function makesnake for
%each condition.
%For example, in the following way:

%% SNAKES
%% UST
folder = 'D:\DatFiles Data\SecondBump_tUni\Data';
files = dir(strcat(folder,'\**\*withall.mat'));
SelArea = 'V1';

% Visual Right
SelConds = [2 3 4];
SelTrialOutcome = [1 0 0];
[frpertrial_all, tempsms, FRz_all, RT_all] = JPpsth_fig6c(files, SelArea, SelConds, SelTrialOutcome);
makesnake(FRz_all,tempsms)
%
SelTrialOutcome = [0 1 0];
[frpertrial_all, tempsms, FRz_all, RT_all] = JPpsth_fig6c(files, SelArea, SelConds, SelTrialOutcome);
makesnake(FRz_all,tempsms)
%
SelTrialOutcome = [0 0 1];
[frpertrial_all, tempsms, FRz_all, RT_all] = JPpsth_fig6c(files, SelArea, SelConds, SelTrialOutcome);
makesnake(FRz_all,tempsms)

% Visual Left
SelConds = [5 6 7];
SelTrialOutcome = [1 0 0];
[frpertrial_all, tempsms, FRz_all, RT_all] = JPpsth_fig6c(files, SelArea, SelConds, SelTrialOutcome);
makesnake(FRz_all,tempsms)
%
SelTrialOutcome = [0 1 0];
[frpertrial_all, tempsms, FRz_all, RT_all] = JPpsth_fig6c(files, SelArea, SelConds, SelTrialOutcome);
makesnake(FRz_all,tempsms)
%
SelTrialOutcome = [0 0 1];
[frpertrial_all, tempsms, FRz_all, RT_all] = JPpsth_fig6c(files, SelArea, SelConds, SelTrialOutcome);
makesnake(FRz_all,tempsms)

% Catch
SelConds = [1];
SelTrialOutcome = [1 1 1];
SelMouseChoice = [1 0 0];
[frpertrial_all, tempsms, FRz_all, RT_all] = JPpsth_fig6c(files, SelArea, SelConds, SelTrialOutcome, SelMouseChoice);
makesnake(FRz_all,tempsms)
%
SelMouseChoice = [0 1 0];
[frpertrial_all, tempsms, FRz_all, RT_all] = JPpsth_fig6c(files, SelArea, SelConds, SelTrialOutcome, SelMouseChoice);
makesnake(FRz_all,tempsms)
%
SelMouseChoice = [0 0 1];
[frpertrial_all, tempsms, FRz_all, RT_all] = JPpsth_fig6c(files, SelArea, SelConds, SelTrialOutcome, SelMouseChoice);
makesnake(FRz_all,tempsms)

%% MST
folder = 'D:\DatFiles Data\SecondBumpTemp';
files = dir(strcat(folder,'\**\*withall.mat'));
SelArea = 'V1';
% Visual Right
SelConds = [2 3 4];
SelTrialOutcome = [1 0 0];
[frpertrial_all, tempsms, FRz_all, RT_all] = JPpsth_fig6c(files, SelArea, SelConds, SelTrialOutcome);
makesnake(FRz_all,tempsms)
%
SelTrialOutcome = [0 1 0];
[frpertrial_all, tempsms, FRz_all, RT_all] = JPpsth_fig6c(files, SelArea, SelConds, SelTrialOutcome);
makesnake(FRz_all,tempsms)
%
SelTrialOutcome = [0 0 1];
[frpertrial_all, tempsms, FRz_all, RT_all] = JPpsth_fig6c(files, SelArea, SelConds, SelTrialOutcome);
makesnake(FRz_all,tempsms)

% Visual Left
SelConds = [5 6 7];
SelTrialOutcome = [1 0 0];
[frpertrial_all, tempsms, FRz_all, RT_all] = JPpsth_fig6c(files, SelArea, SelConds, SelTrialOutcome);
makesnake(FRz_all,tempsms)
%
SelTrialOutcome = [0 1 0];
[frpertrial_all, tempsms, FRz_all, RT_all] = JPpsth_fig6c(files, SelArea, SelConds, SelTrialOutcome);
makesnake(FRz_all,tempsms)
%
SelTrialOutcome = [0 0 1];
[frpertrial_all, tempsms, FRz_all, RT_all] = JPpsth_fig6c(files, SelArea, SelConds, SelTrialOutcome);
makesnake(FRz_all,tempsms)

% Catch
SelConds = [1];
SelTrialOutcome = [1 1 1];
SelMouseChoice = [1 0 0];
[frpertrial_all, tempsms, FRz_all, RT_all] = JPpsth_fig6c(files, SelArea, SelConds, SelTrialOutcome, SelMouseChoice);
makesnake(FRz_all,tempsms)
%
SelMouseChoice = [0 1 0];
[frpertrial_all, tempsms, FRz_all, RT_all] = JPpsth_fig6c(files, SelArea, SelConds, SelTrialOutcome, SelMouseChoice);
makesnake(FRz_all,tempsms)
%
SelMouseChoice = [0 0 1];
[frpertrial_all, tempsms, FRz_all, RT_all] = JPpsth_fig6c(files, SelArea, SelConds, SelTrialOutcome, SelMouseChoice);
makesnake(FRz_all,tempsms)