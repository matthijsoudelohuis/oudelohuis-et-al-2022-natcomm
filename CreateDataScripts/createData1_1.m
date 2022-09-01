%%
params.Experiments         = {'ChangeDetectionConflictDecor' 'VisOnlyPsychophysics' 'BehaviorPsychophysics'};
params.exAnimals           = {'1009' '2034' '2003'};
params.exSessions          = {'2019-03-08_14-39-42' '2019-11-16_15-28-00' '2018-01-17_13-52-00'};

%% Get data:
[Data]                  = MOL_GetData('E:','CHDET',params.Experiments,[],params.exSessions,{'sessionData' 'trialData'});
sessionData             = Data.sessionData;
trialData               = Data.trialData;
trialData               = MOL_RemoveLastnTrials(trialData,20); %Remove last 20 trials:

%% Save the dataset:
save('E:\Matlab\MOL_Analysis\AOudeLohuisetal_2022_NatComms\1Behavior\Data1_1.mat','sessionData','trialData')
fprintf('Dataset: %d sessions, %d trials\n',length(sessionData.session_ID),length(trialData.session_ID));
