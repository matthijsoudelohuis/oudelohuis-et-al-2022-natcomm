% Get arguments:
[Data] = MOL_GetData('E:','CHDET',{'OptoOnly'},{'2009'},{'2018-08-24_11-09-42'},{'sessionData' 'trialData' 'spikeData'});
sessionData         = Data.sessionData;
trialData           = Data.trialData;
spikeData           = Data.spikeData;

TimeWindow          = [sessionData.t_start sessionData.t_stop];
selectedchannel         = 27+32;

%% Get the data:
RawDataDir          = 'K:\Data\CHDET\RawData\2009\2018-08-24_11-09-42\RawData';
%Parameters for Butterworth filter
prm.UseButter               = 0;
%Parameters for Kaiser filter
prm.UseKaiser               = 0;
%Parameters for resampling
prm.UseResampling           = 0;    %Resample data to lower frequency

%Get the data:
lfpData                     = MOL_extractLFP(RawDataDir,selectedchannel,TimeWindow,prm);

%Apply high pass filter:
prm.hp_butter               = 500;  %Low pass filter (Hz) (Butter)
prm.ord_butter              = 4;    %Butterworth filter order
[B_butt,A_butt]             = butter(prm.ord_butter,prm.hp_butter/(lfpData.fs/2),'high');

lfpData.hpsignal{1}         = filtfilt(B_butt,A_butt,lfpData.signal{1});
 
save('E:\Matlab\MOL_Analysis\AOudeLohuisetal_2022_NatComms\5Opto\Data5_1.mat')