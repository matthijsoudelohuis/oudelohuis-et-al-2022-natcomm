params.Experiments          = {'ChangeDetectionConflictDecor' 'VisOnlyPsychophysics' {'BehaviorConflict' 'BehaviorPsychophysics'}};

%% Get the data:
sessionData = struct(); trialData = struct();
for iExp = 1:length(params.Experiments)
    [Data]                  = MOL_GetData('E:','CHDET',params.Experiments{iExp},{},{},{'sessionData' 'trialData'});
    sessionData             = AppendStruct(sessionData,Data.sessionData);
    trialData               = AppendStruct(trialData,Data.trialData);
end
trialData               = MOL_RemoveLastnTrials(trialData,20); %Remove last 20 trials

%% Save the dataset:
save('E:\Matlab\MOL_Analysis\AOudeLohuisetal_2022_NatComms\1Behavior\Data1_2.mat','sessionData','trialData')
fprintf('Dataset: %d sessions, %d trials\n',length(sessionData.session_ID),length(trialData.session_ID));
