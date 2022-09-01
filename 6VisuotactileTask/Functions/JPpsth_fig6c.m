function  [frpertrial_all, tempsms, FRz_all, RT_all, findem, rawspikes] = JPpsth_fig6c(files, SelArea, SelConds, SelTrialOutcome, varargin)
%07-10-2020: added varargin

%Outputs frpertrial (KxTxN), tempsms 
%Also outputs:
% FRz_all : mean across trials zscored per neuron: NxT, only with the neurons used (!don't use this index, for that use output "findem").
% RT_all : reaction times of the trials used. Matches K in frpertrial.
% findem : transforms FRz_all indexes into real indexes: iClu, iSess.
% Average PSTH graph if asked.
% Snake plots graph if asked.
%Inputs:
%   - files : file directory (structure with name, folder, etc.). Example:
%       folder = 'D:\DatFiles Data\SecondBump_tUni\';
%       files = dir(strcat(folder,'\**\*withall.mat'));
%       It needs to direct to files with trialData, events_ts, events_ttl,
%       clus and spikes_clu_ts.
%   - SelArea : select area. V1, PPC, RL, ALL.
%   - SelConds : conditions to utilize. Ex : [2 3]
%   - SelTrialOutcome : boolean vector indicating to include : 
%       [correctResponse firstIncorrect noResponse] Ex: [1 0 0]

% 07/10/2020 : There are many lines that can be pruned out. Can be made
% less messy. Still work to do but will that be useful?

% JPpsth function
    %loop through sessions
    %hardwired input for selection (think about grouping them into a structure)
    %manual variable input for selection (add to structure)
    %function psthselect outputs a filter (input trialData, events, selection criteria)
    %uses filter to build psth
    
    %Ideally 3 filters, each one needs a SelStruct and builds subfilters:
    % cluster filter : subselects neurons (area, fsu, 2b, handpicked etc.)
    % condition filter : subselects event conditions
    % trial filter : subselects trials (multisensory, outcome, side, etc.)
    
    extran = length(varargin);
    switch extran
        case 0
            SelMouseChoice = [1 1 1]; 
        case 1
            SelMouseChoice = varargin{1};
        otherwise
            error('Unexpected inputs')
    end

    
    kk=0;
    FRz_all = [];
    RT_all = [];
    frpertrial_all = cell(1,length(files));
    rawspikes = cell(1,length(files));
    win = [-1000 1200]; %for pietro data: -400 800
    Removepre240 = false;
    

    
for iSess=1:length(files)
    clearvars -except iSess files kk win FRz_all sortmi* RT* frpertrial_all SelArea SelConds SelTrialOutcome SelMouseChoice findem Removepre240 rawspikes
    load(fullfile(files(iSess).folder,files(iSess).name),'spikes_clu_ts','events_ts','events_ttl','trialData','clus') 
    tempsms = win(1):win(2);
    
      % Bin?
        BinSize         = 10; %in ms
        win = [win(1) win(2)-mod(length(win(1):win(2)),BinSize)]; 
      %Other settings :
        boolHalf        = true;
        FilterSize      = 50; %standard deviation of gaussian in ms
  zero_index = find(tempsms>=0, 1, 'first');        

%Select clusters   
 filtclus1 = ones(1,length(spikes_clu_ts)); filtclus2=filtclus1; filtclus3=filtclus1; filtclus4=filtclus1;

    if isfield(clus,'grade')
        %filtclus1 = strcmp(clus.grade,'3') | strcmp(clus.grade,'2');
    end
    if isfield(clus,'p2tcut')
        filtclus2 = ~(clus.p2tcut);
    end
    if isfield(clus,'percviol') && isfield(clus,'fsursu')
        viol1 = (clus.percviol<1.1 & clus.fsursu~=2); %Refr violations <=1.0% for RSU or unknown
        viol2 = (clus.percviol<1.6 & clus.fsursu==2); %Refr violations <=1.5% for FSU 
        filtclus3 = viol1|viol2;
    end
    if isfield(clus,'falsepositiverate')
        %filtclus4 = (clus.falsepositiverate<=2);
    end
%         if isfield(clus,'bump1')
%           filtclus4 =  (clus.fsursu==2); %(clus.bump1==1 | clus.bump2==1); %(clus.fsursu~=1)
%         end

    clusel = 1:length(spikes_clu_ts); 
    clusel = clusel(strcmp(clus.area,SelArea) & filtclus1 & filtclus2 & filtclus3 & filtclus4);
    
%     end
    %ALL wont work, add it.
    
    %Select trigger events to use.
    %Trial-wise:
        %TrialType (Modality)
        trialType_catch     = 1;
        trialType_tactile   = 1;
        trialType_visual    = 1;
        trialType_multi     = 1;

        %Stimulus Side (L vs R)
        %ONLY USE FOR CONGRUENT SESSIONS (for now)
        stimSide_L          = 1;
        stimSide_R          = 1;
        stimSide_N          = 1;

        %Mouse choice resp.
        firstRespLFR_L      = SelMouseChoice(1); 
        firstRespLFR_R      = SelMouseChoice(2); 
        firstRespLFR_none   = SelMouseChoice(3);

        %Trial Outcome
        %   Became input 
        correctResponse     = SelTrialOutcome(1);
        firstIncorrect      = SelTrialOutcome(2);
        noResponse          = SelTrialOutcome(3);

        %Optogenetics
        laserON             = 1;
        laserOFF            = 1;

        %Other event sub-select (not stimuli):
        LFR_End             = 0;
        lickLeft            = 0;
        lickRight           = 0;
        %Epoch sub-select:
        only_StimEpoch      = 0;
        only_LFR            = 0;
        only_ITI            = 0;
    %HISTORY
        %Outcome
        boolUsePreviousOutcome  = 0; %THIS FILTER IS OPTIONAL  
            pCorrect            = 0; %only for previous and actual non-catch trials, otherwise mouse doesnt know the outcome
            pError              = 0; %only for i and i-1 non-catch trials, otherwise mouse doesnt know the outcome
            pMiss               = 0; %only for i and i-1 non-catch trials, for comparison.

% Reparse beh 
                LickSidesInLFR   = cell(1,1);
                            trialData.firstRespLFR = nan(length(trialData.stimStart),1);
                            for k = 1:length(trialData.stimStart)
                                    LickSidesInLFR{k} = trialData.lickSide{k}(trialData.lickTime{k}>=trialData.respwinStart(k) & trialData.lickTime{k}<=trialData.respwinEnd(k));
                                    if ~isempty(LickSidesInLFR{k})
                                    trialData.firstRespLFR(k,1) = LickSidesInLFR{k}(1,1);
                                    end
                            end
                            clear LickSidesInLFR
                            trialData.correctResponse = (trialData.firstRespLFR=='L' & trialData.leftCorrect==1)|...
                                                        (trialData.firstRespLFR=='R' & trialData.rightCorrect==1)|...
                                                        (isnan(trialData.firstRespLFR) & trialData.leftCorrect==0 & trialData.rightCorrect==0);
                            trialData.noResponse = (trialData.firstRespLFR~='L' & trialData.firstRespLFR~='R');
                            trialData.firstIncorrect = (trialData.firstRespLFR=='R' & trialData.leftCorrect==1)|...
                                                        (trialData.firstRespLFR=='L' & trialData.rightCorrect==1)|...
                                                        (~isnan(trialData.firstRespLFR) & trialData.leftCorrect==0 & trialData.rightCorrect==0);;

                numtr = length(trialData.rewardTime);   
                tfields = fieldnames(trialData);
                for iField=1:length(tfields)
                    if (isnumeric( trialData.(tfields{iField}) ) || iscell( trialData.(tfields{iField}) )) && length(trialData.(tfields{iField}))>numtr
                            trialData.(tfields{iField}) = trialData.(tfields{iField})(1:numtr);
                        end
                end
                %Compute RT
                trialData.RT = nan(1,length(trialData.stimStart));
                for k = 1:length(trialData.stimStart)
                    LTIW = trialData.lickTime{k}(trialData.lickTime{k}>=trialData.respwinStart(k) & trialData.lickTime{k}<=trialData.respwinEnd(k));
                            if ~isempty(LTIW)
                                trialData.RT(k) = LTIW(1) - trialData.stimStart(k);
                            end
                end

%
% PSTHLORD_SELECT 2 PART
        EV_TS   = events_ts;  
                idx_stims = find(ismember(events_ttl, [1:49 99])==1); %indexes of stims. can do with str too
                %Subselect trials
                %modality
                Filt1 = zeros(1,length(EV_TS));
                if trialType_catch
                  SCAT=zeros(1,length(EV_TS));
                  SCAT(idx_stims(trialData.trialType=='C'))=1;
                  Filt1 = Filt1|SCAT;
                  %Filt1 = Filt1|strncmp(events_str,'catch',5)  ;
                end
                if trialType_tactile
                  STAC=zeros(1,length(EV_TS));
                  STAC(idx_stims(trialData.trialType=='T'))=1;
                  Filt1 = Filt1|STAC;
                  %Filt1 = Filt1|strncmp(events_str,'tactile',6) ;
                end
                if trialType_visual
                  SVIS=zeros(1,length(EV_TS));
                  SVIS(idx_stims(trialData.trialType=='V'))=1;
                  Filt1 = Filt1|SVIS;
                %   Filt1 = Filt1|strncmp(events_str,'visual',6) ;
                end
                if trialType_multi
                  SMUL=zeros(1,length(EV_TS));
                  SMUL(idx_stims(trialData.trialType=='M'))=1;
                  Filt1 = Filt1|SMUL;
                  %Filt1 = Filt1|strncmp(events_str,'multi',5) ;
                end
                %stimulus side
                Filt2 = zeros(1,length(EV_TS));
                if ~isfield(trialData,'stimSide') %for the rec2 version
                    trialData.stimSide = trialData.stimTACSide;
                    %only for congruency, otherwise visual overwrites.
                    trialData.stimSide(strcmp(trialData.stimVISSide,'R'))={'R'};
                    trialData.stimSide(strcmp(trialData.stimVISSide,'L'))={'L'};
                end
                    if stimSide_L
                        SL = zeros(1,length(EV_TS));
                        SL(idx_stims(strcmp(trialData.stimSide,'L')))=1; %as the rank should be the same.
                        Filt2 = Filt2|SL;
                    end
                    if stimSide_R
                        SR = zeros(1,length(EV_TS));
                        SR(idx_stims(strcmp(trialData.stimSide,'R')))=1;
                        Filt2 = Filt2|SR;
                    end
                    if stimSide_L
                        SN = zeros(1,length(EV_TS));
                        SN(idx_stims(strcmp(trialData.stimSide,'N')))=1;
                        Filt2 = Filt2|SN;
                    end
                %Mouse response
                Filt3 = zeros(1,length(EV_TS));
                if firstRespLFR_L
                    RL = zeros(1,length(EV_TS));
                    RL(idx_stims(trialData.firstRespLFR=='L'))=1;
                    Filt3 = Filt3|RL;
                end
                if firstRespLFR_R
                    RR = zeros(1,length(EV_TS));
                    RR(idx_stims(trialData.firstRespLFR=='R'))=1;
                    Filt3 = Filt3|RR;
                end
                if firstRespLFR_none
                    RN = zeros(1,length(EV_TS));
                    RN(idx_stims(isnan(trialData.firstRespLFR)))=1;
                    Filt3 = Filt3|RN;
                end
                %Trial Outcome
                Filt4 = zeros(1,length(EV_TS));
                if correctResponse
                    CR = zeros(1,length(EV_TS));
                    CR(idx_stims(trialData.correctResponse==1))=1;
                    Filt4 = Filt4|CR;
                end
                if firstIncorrect    
                    FI = zeros(1,length(EV_TS));
                    FI(idx_stims(trialData.firstIncorrect==1))=1;
                    Filt4 = Filt4|FI;
                end
                if noResponse        
                    NR = zeros(1,length(EV_TS));
                    NR(idx_stims(trialData.noResponse==1))=1;
                    Filt4 = Filt4|NR;
                end
                %Optogenetics
                Filt5 = zeros(1,length(EV_TS));
                if laserON     
                    LON = zeros(1,length(EV_TS));
                    LON(idx_stims(trialData.laserON==1))=1;
                    Filt5 = Filt5|LON;    
                end
                if laserOFF       
                    LOFF = zeros(1,length(EV_TS));
                    LOFF(idx_stims(trialData.laserON==0))=1;
                    Filt5 = Filt5|LOFF; 
                end
                %Condition-number
                    Filt6 = zeros(1,length(EV_TS));
                    FCON = zeros(1,length(EV_TS));
                    FCON(idx_stims(ismember(trialData.conditionnumber,SelConds)))=1;
                    Filt6 = Filt6|FCON;  

                %HISTORY
                %Previous Outcome
                Filt7 = ones(1,length(EV_TS));
                if boolUsePreviousOutcome
                    Filt7 = zeros(1,length(EV_TS));
                    if pCorrect
                        PCOR = zeros(1,length(EV_TS));
                        aaa=find((trialData.trialType~='C'))';
                        bbb=find((trialData.correctResponse==1)&(trialData.trialType~='C'))';
                        ccc=find([ismember(aaa,bbb)]==1);
                        PCOR(idx_stims(aaa(ccc(1:end-1)+1))) = 1;
                        Filt7 = Filt7|PCOR;
                    end
                    if pError
                        PERR = zeros(1,length(EV_TS));
                        aaa=find((trialData.trialType~='C'))';
                        bbb=find((trialData.firstIncorrect==1)&(trialData.trialType~='C'))';
                        ccc=find([ismember(aaa,bbb)]==1);
                        PERR(idx_stims(aaa(ccc(1:end-1)+1))) = 1;
                        Filt7 = Filt7|PERR;    
                    end
                    if pMiss
                        PMIS = zeros(1,length(EV_TS));
                        aaa=find((trialData.trialType~='C'))';
                        bbb=find((trialData.noResponse==1)&(trialData.trialType~='C'))';
                        ccc=find([ismember(aaa,bbb)]==1);
                        PMIS(idx_stims(aaa(ccc(1:end-1)+1))) = 1;
                        Filt7 = Filt7|PMIS;
                    end
                end

                %%% FOR OPTO STUFF
                %Remove all trials where response was before 240ms.
                %Removing only delayed ones would introduce bias.
                if Removepre240 
                filtur = ~(trialData.RT<0.24);
                Filt8 = zeros(1,length(EV_TS));
                Filt8(idx_stims(filtur))=1;
                else
                Filt8 = ones(1,length(EV_TS));
                end
                %%%
                %%%

                %Merge filters
                Filt = Filt1&Filt2&Filt3&Filt4&Filt5&Filt6&Filt7&Filt8;
                STIMS_TS = EV_TS(Filt);
                %Subselect other event or just stims ?(
                if ~(LFR_End||lickLeft||lickRight) 
                    %By default just use stims
                    EV_TS = STIMS_TS;

                else
                    %Use other events like licks or LFRend in the subselected trials.
                    if lickLeft
                    %LL_TS = EV_TS(find(strcmp(EV_STR,'lickLeft')));
                    LL_TS = EV_TS(events_ttl==16384 | events_ttl==16386);
                    %Only for trying first lick in bout
                    %DTS = diff([0 LL_TS]);
                    %LL_TS = LL_TS(DTS>0.25*1E6);
                    %
                    EV_TS = LL_TS;
                    end
                    if lickRight
                %     LR_TS = EV_TS(find(strcmp(EV_STR,'lickRight')));
                      LR_TS = EV_TS(events_ttl==32768);
                    %Only for trying first lick in bout
                    %DTS = diff([0 LR_TS]);
                    %LR_TS = LR_TS(DTS>0.25*1E6);
                    %
                    EV_TS = LR_TS;
                    end
                end

                if (only_StimEpoch || only_LFR || only_ITI)

                EE = [];
                %Subselect epoch
                if only_StimEpoch
                    SEEV = [];
                    for k=1:length(STIMS_TS)
                        SEEV =  [SEEV EV_TS(EV_TS>(STIMS_TS(k)) & EV_TS<(STIMS_TS(k)+1E6))];
                    end
                    EE = [EE SEEV];
                end

                if only_LFR
                    LFREV = []; 
                    SelTrials = find(ismember(events_ts(idx_stims), STIMS_TS));
                    LFRStarts = (trialData.respwinStart - trialData.stimStart)*1E6 -0.05*1E6;
                    LFREnds = (trialData.respwinEnd - trialData.stimStart)*1E6 +0.05*1E6;    
                    %We give a tolerance of ??ms after LFR for licks triggered by
                    %inside-LFR stimuli
                    for k=1:length(STIMS_TS)
                        JJJ = EV_TS(EV_TS>(STIMS_TS(k)+LFRStarts(SelTrials(k))) & EV_TS<(STIMS_TS(k)+LFREnds(SelTrials(k))));
                        %Onlt first lick
                        if ~isempty(JJJ)
                            JJJ = JJJ(1);
                        end
                        LFREV = [LFREV JJJ]; 
                        clear JJJ
                    end
                    EE = [EE LFREV];
                end

                if only_ITI
                    ITIEV = [];
                    for k=1:length(STIMS_TS)
                        ITIEV = [ITIEV EV_TS(EV_TS>(STIMS_TS(k)-2.5*1E6) & EV_TS<(STIMS_TS(k)-10))]; 
                    end
                    EE = [EE ITIEV];
                end

                EV_TS = unique(EE);
                end

                %Subselect



                %Also have a trialData version for Filt
                if ~(LFR_End||lickLeft||lickRight) 
                    trialDataFilt = ismember(events_ts(idx_stims),EV_TS);
                    trialDataF = trialData;
                    for iField=1:length(tfields)
                        trialDataF.(tfields{iField}) = trialData.(tfields{iField})(trialDataFilt);
                    end
                    
                end
%               
% Make PSTH per neuron 
tempsms_bins = win(1):BinSize:win(2);
frpertrial_all{iSess} = nan(length(EV_TS),length(tempsms_bins),length(spikes_clu_ts));
temps_1ms = win(1):win(2); %for raster
rawspikes{iSess} = nan(length(EV_TS),length(temps_1ms),length(clusel)); %or temps_1ms
clucounter = 0 ;

    for clu = clusel
        clucounter = clucounter + 1;
        STS     = spikes_clu_ts; 
        STS = STS(clu);
        
        %init
        tempsms = win(1):win(2);
        tempsms_bins = win(1):BinSize:win(2);
        frpertrial = nan(length(EV_TS),length(tempsms));
        frpertrial_bins = nan(length(EV_TS),length(tempsms_bins));
        
        %If less than 2 trials it's useless
        if length(EV_TS)<=0  %2 
            FRmean = nan(size(tempsms_bins));
            FRz = nan(size(tempsms_bins));   
            RT=nan(1,length(EV_TS));
            frpertrial = frpertrial_bins;
            tempsms = tempsms_bins;
        else
        %Initialize PSTH stuff
        rast=[];
        rastcon = [];
        raw_psth={};
        Gauss_width = max([11 6*FilterSize+1]);
        kernel      = normpdf(-floor(Gauss_width/2):floor(Gauss_width/2),0,FilterSize);
        if boolHalf
            kernel = pdf('HalfNormal',-floor(Gauss_width/2):floor(Gauss_width/2),0,FilterSize);
        end
        kernel_bins = kernel(1:10:end);
        RT=nan(1,length(EV_TS));
        % Start real PSTH
            for k=1:length(EV_TS) %number of trials
                if ~(LFR_End||lickLeft||lickRight) 
                LTIW = trialDataF.lickTime{k}(trialDataF.lickTime{k}>=trialDataF.respwinStart(k) & trialDataF.lickTime{k}<=(trialDataF.respwinEnd(k)+0.01));
                if ~isempty(LTIW)
                    RT(k) = LTIW(1) - trialDataF.stimStart(k);
                end
                end
                app = find(STS{1,1}>=(EV_TS(k)+win(1)*1E3) & STS{1,1}<=(EV_TS(k)+win(2)*1E3)); %find spikes ts in window
                raw_psth{k} = STS{1,1}(app)-EV_TS(k); %make times relative to event
                t = raw_psth{k}; 
                rast = [rast, t]; %Stores ALL relative spiketimes for PSTH
                %make it in timebase 1ms then convolve w gaussian.
                %frpertrial(k,:) = conv( histc(raw_psth{k}/1E3, tempsms) , kernel, 'same' ) / (size(STS,2) * 1E-3); %can use valid for only valid part of conv    
                %COMB = conv( histc(raw_psth{k}/1E3, tempsms) , kernel, 'valid' );
                rawspikes{iSess}(k,:,clucounter) = histc(raw_psth{k}/1E3, temps_1ms); %or temps_1ms
                COMB_bins = conv( histc(raw_psth{k}/1E3, tempsms_bins) , kernel_bins, 'valid' );
                %frpertrial(k,:) = [nan(1,floor((length(tempsms)-length(COMB))/2)) COMB nan(1,ceil((length(tempsms)-length(COMB))/2))] / (size(STS,2) * 1E-3); %can use valid for only valid part of conv
                frpertrial_bins(k,:) =  [nan(1,floor((length(tempsms_bins)-length(COMB_bins))/2)) COMB_bins nan(1,ceil((length(tempsms_bins)-length(COMB_bins))/2))] / (size(STS,2) * 1E-3);
                %frpertrial_bins(k,:) =  COMB_bins / (size(STS,2) * 1E-3);
% % %                 %%% for opto figure : nans after LFR end
                if Removepre240
                lfrend_t_bins = find(tempsms_bins>=1000*(trialDataF.respwinEnd(k) - trialDataF.stimStart(k)),1);
                if ~isempty(lfrend_t_bins)
                    frpertrial_bins(k,:) = [frpertrial_bins(k,1:lfrend_t_bins) nan(1,length(tempsms_bins)-lfrend_t_bins)];
                end
                end
% % %                 %%% for opto figure
            end
        %     zero_index = find(tempsms>=0, 1, 'first');
        %     frpertrialz =  nanmean(frpertrial(:,1:zero_index),2)  ;

        %For now just pretend it's the usual frpertrial with 1ms bins
        frpertrial = frpertrial_bins;
        tempsms = tempsms_bins;

        FRmean = nanmean(frpertrial,1);
        
%         FRz = (FRmean - nanmean(FRmean))./nanstd(FRmean);
        if isfield(clus,'z_mutri')
%             FRz = (FRmean - clus.z_mutri(clu))./ clus.z_sigtri(clu) ; 
            FRz = (FRmean - clus.z_mu_mol(clu))./ clus.z_sig_mol(clu) ; 
        else
%             error(strcat('No z-score mol in cluster:',num2str(clu),', session:',num2str(iSess)))
%             zero_index = find(tempsms>=0, 1, 'first');
%             FRbl_smooth = nanmean(FRmean(1:zero_index));
%             FR_std = nanstd(FRmean);
%             inin = find(tempsms>=-500,1);
%             A = frpertrial_all{1,1}(:,1:inin,:);
%             B = squeeze(reshape(A,1,size(A,1)*size(A,2),size(A,3)));
%             FR_std = nanstd(B);
%             FRz = (FRmean - FRbl_smooth)./FR_std; %
              FRz = (FRmean - nanmean(FRmean))./nanstd(FRmean);
        end
% % % FRz = FRmean./FRbl_smooth ;

        % Use all trial for z-scoring so that std doesn't equal 0. 
        % The only case then it will equal 0 would be no spikes, in that case we
        % can replace the Inf in the z-score for a nan
        FRz(FRz==inf) = nan;


        end
    
        kk=kk+1;
        FRz_all(kk,:) = FRz;
        findem(kk,1) = clu;
        findem(kk,2) = iSess;
          if exist('frpertrial')
        frpertrial_all{iSess}(:,:,clu) = frpertrial;
          end
          
    end
        
        RT_all{iSess} = RT;
end


end