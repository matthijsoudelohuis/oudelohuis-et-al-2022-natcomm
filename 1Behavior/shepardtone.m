function [tone,Fs] = shepardtone(freq,tones,varargin)
%SHEPARDTONE Calculate Shepard tone for the shepard scale
%	[TONE,FS] = SHEPARDTONE(FREQ,TONES) calculates tone with 
%   center frequency FREQ and numer of tones TONES. The tones are weighted 
%   with a standard Gauss window. TONES, the number of partial tones, 
%   must be odd, otherwise the next lower odd integer is taken. 
%   The resulting tone is stored in TONE.
% 
%   [TONE,FS] = SHEPARDTONE(FREQ,TONES,'PARAM1',val1, 'PARAM2',val2,...)
%   specifies one or more of the following name/value pairs:
%
%     'length'       Specifies the duration of the shepard tone
%                    (default 0.2 sec)
%
%     'samplerate'   Desired samplerate
%                    (default 44100 Hz)
%
%     'weight'       Specifies special weight vector
%                    (default is standard Gauss window)
%
%
% Example 1:
%   Shepard tone of 440 Hz center frequency with 5 partial tones standard
%   Gauss-weighted.
%   [tone,fs] = shepardTone(440,5);
%   sound(tone,fs);
%
% Example 2:
%   A shepard tone is created with sample rate 22100 samples/s and a 
%   special weight vector. The resultig tone lasts for 1 sec.
%   [tone,fs] = shepardTone(440,5,'samplerate',22100,'length',1,...
%              'weight',gausswin(5,1.5));
%   sound(tone,fs);
%
% Reference:
%   Shepard, R. (1964) "Circularity in Judgements of Relative Pitch". 
%   Journal of the Acoustical Society of America 36 (12): 2346-53. 
%
% Frederik Nagel
% Institute of Music Physiology and Musicians' Medicine
% Hanover University of Music and Drama 
% Hanover
% Germany
%
% e-mail: frederik.nagel@hmt-hannover.de
% homepage: http://www.immm.hmt-hannover.de
%
% Oct 6, 2006.

% Odd?
if(mod(tones,2)==0)
    tones = tones-1;
    disp(['Note! Number of tones is changed to ' num2str(tones) '!']);
    
end

%defaults
len = .2;
Fs = 44100;
weight = gausswin(tones);
par = 'x';

if(nargin > 2)
    for i=1:2:size(varargin,2)
        par = lower(varargin{i});
        if ~isnumeric(varargin{i+1})
            error('Wrong parameters!');    
        else
            val = varargin{i+1};
        end
        switch(par)
            case 'length'
                len = val;
            case 'samplerate'
                Fs = val;
            case 'weight'
                weight = val;
                if(length(weight)~=tones)
                    error('Vector of weights has different length than tones.');
                end
        end
    end
end


% Sample rate
T = 1/Fs;
%t = T:T:len;
t = [0:T:len-T];

tones = tones+1;
i=(1:tones-1)'; %All partial tones
% shepard tone
TONE = sum(sin(freq * 2.^(i-tones/2) *(2*pi)*t).* repmat(weight(i),1,length(t)),1);

% Normalize
%tone = TONE/max(TONE);
tone = TONE;