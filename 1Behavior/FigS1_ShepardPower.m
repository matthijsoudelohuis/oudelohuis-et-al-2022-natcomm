%% Shepard tone spectrogram
% This script shows the spectrogram for the auditory stimuli
% Oude Lohuis et al. 2022 Nat Comms

Par.strAuStimType               = 'shepardtone';
Par.intSamplingRate             = 192000;                    %Preferred sampling rate
Par.FreqSpacing                 = 1/256;                    %Spacing between consecutive frequencies in octaves
Par.allOctaves                  = 13:Par.FreqSpacing:14-Par.FreqSpacing; %Init all octaves
Par.vecLoadFreqs                = 2.^Par.allOctaves;        %Init all frequencies
Par.vecLoadFreqs                = round(Par.vecLoadFreqs);  %Round them off to closest integer
% Par.vecLoadFreqs                = 8000:10:16000-10;  %Round them off to closest integer

Par.nFreqs                      = length(Par.vecLoadFreqs); %store number of loaded frequencies for module circularity
Par.ShepardTones                = 5;                        %Number of shepard tones (partial tones above and below)
Par.ShepardWeights              = gausswin(length(Par.vecLoadFreqs) * Par.ShepardTones); %Generate partial tone weights divide for sound range 0-1
Par.allHarmonics                = transpose((13-2):Par.FreqSpacing:(14+2)-Par.FreqSpacing);
Par.vecFreqChange               = 0.5;                      %Vector with possible delta frequency in octave
Par.minSampleDur                = 0.02;                     %Minimum duration of loaded tone that is repeated (too low and the schedule might get empty)
Par.maxSampleDur                = 0.1;                      %Maximum duration of loaded tone that is repeated (this determines resolution for changing ad hoc)

intThisTrial                    = 1;
Par.Stim.vecAuStimInt(intThisTrial) = 1;

%% Figure of weights in compound Shepard tone
exampletones = [64 196];

colors = {[245 145 24] [227 1 245] [120 245 24] [12 170 245]}; colors = cellfun(@(x) x/256,colors,'UniformOutput',false);

figure; set(gcf,'units','normalized','Position',[0.6 0.6 0.22 0.23],'color','w'); hold all;

% patch([13 14 14 13],[0 0 1 1],[0.8 0.8 0.8],'EdgeColor',[0.8 0.8 0.8]);
plot(Par.allHarmonics,Par.ShepardWeights,'k','LineWidth',2)
% bar(Par.allHarmonics,Par.ShepardWeights,'k','LineWidth',0.02)

% idx = 1:256:length(Par.allHarmonics);
% temp = zeros(size(Par.ShepardWeights)); temp(idx) = Par.ShepardWeights(idx);
% plot(Par.allHarmonics,temp,'-','Color',[0.8 0.8 0.8])

idx = exampletones(1):256:length(Par.allHarmonics);
% plot(Par.allHarmonics(idx),Par.ShepardWeights(idx),'.','Color',[0 0.6 0.3],'MarkerSize',25)
plot(Par.allHarmonics(idx),Par.ShepardWeights(idx),'.','Color',colors{2},'MarkerSize',25)
temp = zeros(size(Par.ShepardWeights)); temp(idx) = Par.ShepardWeights(idx);
plot(Par.allHarmonics,temp,'-','Color',colors{2})
idx = exampletones(2):256:length(Par.allHarmonics);
plot(Par.allHarmonics(idx),Par.ShepardWeights(idx),'.','Color',colors{3},'MarkerSize',25)
temp = zeros(size(Par.ShepardWeights)); temp(idx) = Par.ShepardWeights(idx);
plot(Par.allHarmonics,temp,'-','Color',colors{3})
xlim([11 16])
ylim([0.001 1])

%% Make small figures for the weights of four cardinal tones
exampletones = [1 64 128 196];

figure; set(gcf,'units','normalized','Position',[0.06 0.6 0.46 0.12],'color','w'); hold all;

for i = 1:4
    subplot(1,4,i); hold all;
    plot(Par.allHarmonics,Par.ShepardWeights,'k','LineWidth',2)
    idx = exampletones(i):256:length(Par.allHarmonics);
    temp = zeros(size(Par.ShepardWeights)); temp(idx) = Par.ShepardWeights(idx);
    plot(Par.allHarmonics(idx),Par.ShepardWeights(idx),'.','Color',colors{i},'MarkerSize',25)
    plot(Par.allHarmonics,temp,'-','Color',colors{i})
    xlim([11 16])
    ylim([0.001 1])
end

%% Make spectrogram of example auditory trial change:
exFreqs                                 = [64 196 196 196-16];

clear varAudio 
varAudio(length(Par.vecLoadFreqs))          = struct();

vecSoundChange          = [];

% for freq = find(idx)
for freq = exFreqs
    varAudio(freq).oct                  = Par.allOctaves(freq);
    varAudio(freq).freq                 = Par.vecLoadFreqs(freq);
    Par.Stim.vecFreq(intThisTrial)      = Par.vecLoadFreqs(freq);
    tempvar                             = load_auditory(Par,intThisTrial);
    stimDur                             = size(tempvar.vecSound,2)/Par.intSamplingRate;
    vecSound                            = repmat(tempvar.vecSound(1,:),1,floor(1/stimDur));
    vecSoundChange                      = [vecSoundChange vecSound];
end

%Compute the short-time Fourier transform
spacing = 1/256;
specOct = 10:spacing:17;
specF = 2.^specOct;
[S,F,T,P] = spectrogram(vecSoundChange,1000,100,specF,192000);

%trim very low powers
% P(P<1.9618e-12) = 1.9618e-12;
P(P<1e-12) = 1e-12;

% Make the figure
figure; set(gcf,'units','normalized','Position',[0.6 0.6 0.43 0.25],'color','w')
% imagesc(10*log10(psdx_all)); hold on;
imagesc(P); hold on;
% hotmap = parula(1200);
% hotmap = getPyPlot_cMap('CMRmap');
% hotmap = getPyPlot_cMap('Reds');
hotmap = [1:-0.01:0; 1:-0.01:0; 1:-0.01:0]'; hotmap(:,1) = hotmap(:,1).^0.25; %hotmap(:,3) = hotmap(:,3).^0.5;
cmap = colormap(hotmap); % caxis([-max(max(psdx_all)) 0])
colormap(cmap.^4);
caxis([1e-12 1e-4])

title('Spectrogram','FontSize',15)
resol = 128;
set(gca,'YTick',1:resol:length(specF), 'YTickLabels',specOct(1:resol:end),'FontSize', 12)
set(gca,'XTick',find(T>1,1), 'XTickLabels','Frequency Change','FontSize', 10)
xlabel('Time (s)')
ylabel('Frequency (Oct)')

c = colorbar;
c.Ticks = [0.001 0.5 1]*1e-4;
c.TickLabels = {'10^{-12}' '10^{-8}' '10^{-4}'};

%% Construct spectrogram of all possible tones:
nSamples                        = 10000;
psdx_all        = NaN(Par.nFreqs,nSamples/2+1);
SPL             = NaN(Par.nFreqs,1);
for iFreq = 1:Par.nFreqs
    
    Par.Stim.vecFreq(intThisTrial)    = Par.vecLoadFreqs(iFreq);
    varaudio            = load_auditory(Par,intThisTrial);
    sound               = repmat(varaudio.vecSound(1,:),1,5000);
    sound               = sound(1:nSamples);
%     sound               = filterA(sound,Par.intSamplingRate); 
%     SPL(iFreq)          = spl(sound,'air');
    x                   = sound + randn(size(sound))/10000; %Add noise to facilitate log10 imagesc
    N                   = length(x);
    xdft                = fft(x);
    xdft                = xdft(1:N/2+1);
    psdx                = (1/(Par.intSamplingRate*N)) * abs(xdft).^2;
    psdx(2:end-1)       = 2*psdx(2:end-1);
%     psdx                = smooth(psdx,10,'sgolay',2)';
    psdx                = smooth(psdx,10)';
    psdx_all(iFreq,:)   = psdx;
end

freq                    = 0:Par.intSamplingRate/nSamples:Par.intSamplingRate/2;

%Make figure:
figure; set(gcf,'units','normalized','Position',[0.2 0.1 0.25 0.28],'color','w')
imagesc(psdx_all); hold on; 

hotmap = parula(1000);
cmap = colormap(hotmap); 
colormap(cmap.^1.2);
caxis([1e-12 1e-4])

[~,lowIdx]              = min(abs(freq-Par.vecLoadFreqs(1)));
[~,highIdx]             = min(abs(freq-Par.vecLoadFreqs(end)));
plot([lowIdx highIdx highIdx lowIdx lowIdx],[Par.nFreqs Par.nFreqs 1 1 Par.nFreqs],'k','LineWidth',2);

title('Power Spectrum','FontSize',12)

resol = 64;

set(gca, 'YTick', 1:64:256, 'YTickLabels',sprintfc('2^{%2.2f}',Par.allOctaves([1:64:256])),'FontSize', 12)

tickpos = [find(freq>=2^12,1) find(freq>=2^13,1) find(freq>=2^14,1) find(freq>=2^15,1) find(freq>=2^16,1)];
tickposlabels = {'4098' '8192' '16384' '32768' '65536'};
set(gca, 'XTick',tickpos, 'XTickLabels',tickposlabels,'FontSize', 10)
xlabel('Frequency (kHz)')
ylabel('Ground Tone (kHz)')

