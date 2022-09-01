params.minPhotostimpower   = 2;    %Only include session in which power was >2mW total

%% Load all data:
%normal trialData (as in behavioral sessions, i.e. no extra trials inserted)
[Data] = MOL_GetData('E:','CHDET',{'ChangeDetectionConflict' 'VisOnlyTwolevels'},{'2019' '2020' '2021' '2022' '2023' '2026' '2027' '2030' '2031' '2028' '2029' '2034' '2035'},[],{'sessionData' 'trialData'});
sessionData     = Data.sessionData;
trialData       = Data.trialData;

%% Remove last 20 trials:
trialData = MOL_RemoveLastnTrials(trialData,20);

%% Filter sessions with optogenetic manipulation in V1:
sesids            = unique(sessionData.session_ID(sessionData.UseOpto==1 & ...
    (strcmp(sessionData.PhotostimArea,'V1') | strcmp(sessionData.PhotostimArea,'S1')) & ...
    sessionData.Photostimpower >= params.minPhotostimpower));
fprintf('Selected %d/%d sessions with optical manipulation in V1 or S1\n',length(sesids),length(sessionData.session_ID));
[sessionData, trialData]        = MOL_getTempPerSes(sesids,sessionData,trialData);

%% Filter out multisensory sessions that have too low visual or auditory performance:
nSessions           = length(sessionData.session_ID);
visperf             = NaN(nSessions,1);
auperf              = NaN(nSessions,1);
for iSes = 1:nSessions
    if ismember(sessionData.Experiment(iSes),{'ChangeDetectionConflict' 'BehaviorConflict'})
        sesid                                                       = sessionData.session_ID(iSes);
        [tempsessionData,temptrialData]                             = MOL_getTempPerSes(sesid,sessionData,trialData);%Get the sessionData for each session individually:
        
        vistrialidx     = strcmp(temptrialData.trialType,'X') & temptrialData.visualOriChangeNorm==3 & temptrialData.hasphotostim~=1;
        visperf(iSes)   = sum(temptrialData.correctResponse(vistrialidx)) / sum(vistrialidx);
        autrialidx      = strcmp(temptrialData.trialType,'Y') & temptrialData.audioFreqChangeNorm==3 & temptrialData.hasphotostim~=1;
        auperf(iSes)    = sum(trialData.correctResponse(autrialidx)) / sum(autrialidx);
    end
end

sesids              = sessionData.session_ID(~(visperf<0.3 | auperf<0.3));
fprintf('Removed %d/%d sessions with low behavioral accuracy\n',length(sessionData.session_ID)-length(sesids),length(sessionData.session_ID));
[sessionData, trialData]        = MOL_getTempPerSes(sesids,sessionData,trialData);

%% save the data:
save('E:\Matlab\MOL_Analysis\AOudeLohuisetal_2022_NatComms\5Opto\Data5_3.mat','sessionData','trialData')

