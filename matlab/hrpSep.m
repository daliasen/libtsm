function [xHarm,xPerc,xRes,sideinfo] = hrpSep(x,parameter,sideinfo)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Name: hrpSep
% Date: 03-2022
% Programmer: daliasen
%
% Seperates a given audio signal into a harmonic, a residual, and a
%   percussive component according to the paper:
%
% Jonathan Driedger, Meinard Müller, Sascha Disch
% Extending Harmonic–Percussive Separation of Audio Signals.
% Proceedings of the International Conference on Music Information
% Retrieval (ISMIR), Taipei, Taiwan, 2014, pp. 611–616.
%
% The code is adapted from function hpSep.
%
% Input:  x                 input signal.
%         parameter.
%          anaHop           the stft hop size of the analysis window.
%          win              the stft analysis window used for windowing the
%                           input signal.
%          zeroPad          number of zeros that should be padded to the
%                           window to increase the fft size and therefore
%                           the frequency resolution.
%          filLenHarm       length of the median filter in time direction.
%          filLenPerc       length of the median filter in frequency
%                           direction.
%          beta             separation factor.
%
% Output: xHarm             the harmonic component of the input signal x.
%         xPerc             the percussive component of the input signal x.
%         xRes              the residual component of the input signal x.
%
%         sideinfo.
%            hpSep.stftAnaHop
%            hpSep.win
%            hpSep.zeroPad
%            hpSep.filLenHarm
%            hpSep.filLenPerc
%            hpSep.beta
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% check parameters
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
if nargin < 3
    sideinfo = [];
end
if nargin < 2
    parameter = [];
end
if nargin < 1
    error('Please specify input data x.');
end

if ~isfield(parameter,'anaHop')
    parameter.anaHop = 256;
end
if ~isfield(parameter,'win')
    parameter.win = win(1024,2); % hann window
end
if ~isfield(parameter,'zeroPad')
    parameter.zeroPad = 0;
end
if ~isfield(parameter,'filLenHarm')
    parameter.filLenHarm = 10;
end
if ~isfield(parameter,'filLenPerc')
    parameter.filLenPerc = 10;
end
if ~isfield(parameter,'beta')
    parameter.beta = 2;
end


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% some pre calculations
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
anaHop = parameter.anaHop;
w = parameter.win;
zeroPad = parameter.zeroPad;
filLenHarm = parameter.filLenHarm;
filLenPerc = parameter.filLenPerc;
beta = parameter.beta;
numOfChan = size(x,2);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% harmonic-percussive separation
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
xHarm = zeros(size(x,1),numOfChan);      % Initialize output
xRes = zeros(size(x,1),numOfChan);      % Initialize output
xPerc = zeros(size(x,1),numOfChan);      % Initialize output

for c = 1 : numOfChan                   % loop over channels
xC = x(:,c);

% stft
parStft.anaHop = anaHop;
parStft.win = w;
parStft.zeroPad = zeroPad;
spec = stft(xC,parStft);
magSpec = abs(spec);

% harmonic-percussive separation
magSpecHarm = medianFilter(magSpec,filLenHarm,2);
magSpecPerc = medianFilter(magSpec,filLenPerc,1);

maskHarm = magSpecHarm > beta * magSpecPerc;
maskPerc = magSpecPerc >= beta * magSpecHarm;
maskRes = 1 - (maskHarm + maskPerc);

specHarm = maskHarm .* spec;
specRes = maskRes .* spec;
specPerc = maskPerc .* spec;

% istft
parIstft.synHop = parameter.anaHop;
parIstft.win = parameter.win;
parIstft.zeroPad = parameter.zeroPad;
parIstft.numOfIter = 1;
parIstft.origSigLen= length(x);
xHarmC = istft(specHarm,parIstft);
xResC = istft(specRes,parIstft);
xPercC = istft(specPerc,parIstft);

xHarm(:,c) = xHarmC;
xRes(:,c) = xResC;
xPerc(:,c) = xPercC;
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% update sideinfo
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
sideinfo.hpSep.stftAnaHop = parameter.anaHop;
sideinfo.hpSep.win = parameter.win;
sideinfo.hpSep.zeroPad = parameter.zeroPad;
sideinfo.hpSep.filLenHarm = parameter.filLenHarm;
sideinfo.hpSep.filLenPerc = parameter.filLenPerc;

end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% median filter
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function Y = medianFilter(X,len,dim)

s = size(X);
Y = zeros(s);

switch dim
    case 1
        XPadded = [zeros(floor(len/2),s(2));X;zeros(ceil(len/2),s(2))];
        for i = 1 : s(1)
            Y(i,:) = median(XPadded(i:i+len-1,:),1);
        end

    case 2
        XPadded = [zeros(s(1),floor(len/2)) X zeros(s(1),ceil(len/2))];
        for i = 1 : s(2)
            Y(:,i) = median(XPadded(:,i:i+len-1),2);
        end

    otherwise
        error('unvalid dim.')
end

end
