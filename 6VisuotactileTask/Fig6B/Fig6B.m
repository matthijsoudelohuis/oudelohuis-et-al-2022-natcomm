% This code generates Figure 6b from oude Lohuis et al., 2022.
% Multinomial logistic regression is used to fit psychometric curves on behavioral data from Task B. 
% For any questions contact Umberto Olcese at u.olcese@uva.nl

%% Fit data per mouse and generate fit files
clear all
%Insert directory where OudeLohuisetal_2022_NatComms is saved:
maindirectory = 'C:\Users\jeanl\surfdrive\Shared\Manuscript - 2nd Bump\Submission - NatComm\rev1';

%% For unisensory (UST) mice 
folder = strcat(maindirectory,'\OudeLohuisetal_2022_NatComms\6VisuotactileTask\Data6_1\UST');
files = dir(strcat(folder,'\**\All\*_all.mat'));

%Generate fit files
for iMuis = 1:length(files)
    cd(fullfile(files(iMuis).folder))
    load(fullfile(files(iMuis).folder,files(iMuis).name)) 
    PieTaskBeh_fittingfig1
    MLogisticReg_fig1
    clearvars -except folder files iMuis
end

%% For multisensory (MST) mice
folder = strcat(maindirectory,'\OudeLohuisetal_2022_NatComms\6VisuotactileTask\Data6_1\MST');
files = dir(strcat(folder,'\**\all\*_all.mat'));
for iMuis = 1:length(files)
    cd(fullfile(files(iMuis).folder))
    load(fullfile(files(iMuis).folder,files(iMuis).name)) 
    MLogisticReg_bothbehs_fig1
    clearvars -except folder files iMuis
end

%% Use fit files to plot fits for all mice
Plot_Fig6B