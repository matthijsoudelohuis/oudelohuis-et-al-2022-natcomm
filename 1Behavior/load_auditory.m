%% Load the auditory stimulus
% Return varAudio, a struct which contains the information about the audio
% stimulus
function varAudio = load_auditory(Par,intThisTrial)

if isfield(Par,'strAuStimType')    % Generate Audio
    
    varAudio = []; % variable for storing auditory stimulus
    if strcmp(Par.strAuStimType,'bandpassfreq')
        Stimdur     =   5; %seconds (Use sufficient time for stimulus, stimulus is stopped by task script)
    elseif strcmp(Par.strAuStimType,'shepardtone')
        %Approach: this part will increase the number of cycles (full oscillations) to get the amount of
        %cycles that matches with the sampling frequency
        %(then no artefact when repeating the stimulus)
        minfreq       = Par.Stim.vecFreq(intThisTrial)*2^-((Par.ShepardTones-1)/2);
        
        cycles        = 1; %start with one oscillation
        Stimdur       = NaN; %Init stimdur
        while mod(Stimdur,1/Par.intSamplingRate)~=0 %keep increasing cycle until mod=0
            cycles = cycles+1;
            Stimdur             = 1/minfreq*cycles;
            
            if Stimdur>Par.maxSampleDur %If the number of cycles 
                maxcycles     = round(Par.maxSampleDur/(1/minfreq));
                [~,cycles]       = min(mod([1:maxcycles]*1/minfreq,1/Par.intSamplingRate));
                break; %terminate from while loop
            end
        end
        Stimdur = 1/minfreq*cycles;
        
        while Stimdur<Par.minSampleDur %keep increasing stim dur until stimdur long enough:
            Stimdur             = Stimdur*2;
        end
    end
    
    %Generate auditory stimulus
    if strcmp(Par.strAuStimType,'bandpassfreq')
        vecSoundOneSide = white_noise_band_ramp(Stimdur, Par.intSamplingRate,...
            Par.Stim.vecFreq(intThisTrial)-Par.FreqBandwidth/2,...
            Par.Stim.vecFreq(intThisTrial)+Par.FreqBandwidth/2,...
            Par.dblAudioRampDur);
    elseif strcmp(Par.strAuStimType,'shepardtone')
        window               = Par.ShepardWeights(find(Par.vecLoadFreqs == Par.Stim.vecFreq(intThisTrial)):length(Par.vecLoadFreqs):Par.ShepardTones*length(Par.vecLoadFreqs));
        window               = window/2.5;
        vecSoundOneSide      = shepardtone(Par.Stim.vecFreq(intThisTrial),Par.ShepardTones,'length',Stimdur,'samplerate',Par.intSamplingRate,'weight',window);
    end
    
    %Convert intensity values (0-1) to output in dB:
    if strcmp(Par.strAuStimType,'bandpassfreq')
        mappingintenstodb = 1/(2^((1-Par.Stim.vecAuStimInt(intThisTrial))*10));
    elseif strcmp(Par.strAuStimType,'shepardtone')
        mappingintenstodb = 1/(3.3^((1-Par.Stim.vecAuStimInt(intThisTrial))*10));
    end
    
    % Set auditory intensity for both sides and store in variable varaudio
    varAudio.vecSound = [vecSoundOneSide * mappingintenstodb; vecSoundOneSide * mappingintenstodb];
    
else
    % retrieve side for current trial
    charSide = Par.Stim.vecSide(intThisTrial);
    varAudio = []; % variable for storing auditory stimulus
    
    %Generate auditory stimulation (white noise)
    vecSoundOneSide = white_noise_band_ramp(Par.dblSecsStimDur, Par.intSamplingRate,...
        Par.vecNoiseRange(1),Par.vecNoiseRange(2),Par.dblAudioRampDur); % in other file
    % right
    if charSide == 'R'
        varAudio.vecSound = [vecSoundOneSide; zeros(1,numel(vecSoundOneSide))];
        % left
    elseif charSide == 'L'
        varAudio.vecSound = [zeros(1,numel(vecSoundOneSide)); vecSoundOneSide];
        % center
    elseif charSide == 'C'
        varAudio.vecSound = [vecSoundOneSide; vecSoundOneSide];
    end
    varAudio.intAudioChannel = Par.intAudioChannel;
    
end


end
