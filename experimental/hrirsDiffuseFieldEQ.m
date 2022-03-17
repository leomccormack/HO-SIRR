function [eqTaps] = hrirsDiffuseFieldEQ(hrirs, MINPHASE, gridWeights)
%hrirsDiffuseFieldEQ Calculate diffuse field (common transfer function) EQ.
%   INPUTS
%       hrirs : [len x 2 x grid]
%       MINPHASE : bool (optional, default true)
%       gridWeights : [grid, 1] (optional, default [] assumes regular grid)
%   OUTPUT
%       eqTaps : Filter taps [len x 1]
%       
% Chris Hold 2021,2022

if nargin < 3
    gridWeights = [];
end
if nargin < 2
    MINPHASE = true;
end

numTaps = size(hrirs, 1);
numGrid = size(hrirs, 3);

% Check if grid is given
if isempty(gridWeights)
    gridWeights = (4*pi) / numGrid;  % Close enough in this case
else
    assert(size(gridWeights, 1) == numGrid)
    assert(size(gridWeights, 2) == 1)
end


% FD transform
nfft = 16*numTaps;  % interpolate
H = fft(hrirs, nfft, 1);
Hs = H(1:nfft/2+1, :, :);


% weighted RMS
Havg = sqrt(sum(reshape(gridWeights, 1, 1, []) .* abs(Hs).^2, 3) / (4*pi));

% Avg (left, right)
Havg = mean(Havg, 2);

% Smoothing
HavgSmooth = Havg;
for bin = 2:nfft/2+1
    if ~mod(bin, 2)
        avgidx = (bin - bin/2 : min(3/2 * bin, nfft/2+1));
    end
    win = hann(length(avgidx));
    % weighted average
    HavgSmooth(bin) = sum(win .* Havg(avgidx)) / sum(win);
end

% frequency mask for lo and hi, avoid inversion there
freqWeight = ones(nfft/2+1, 1);
wLo = hann(nfft/64+1);
freqWeight(1:nfft/128 + 1) = wLo(1:nfft/128 + 1);
wHi = hann(nfft/2 + 1);
freqWeight(end-nfft/4:end) = wHi(end-nfft/4:end);

% Avoid excessive filters
HavgSmooth = HavgSmooth / mean(HavgSmooth);

% frequency weighted inversion 
HinvWeighted = (freqWeight).*(1 ./ HavgSmooth) + ...
               (1-freqWeight) .* ones(size(HavgSmooth));



% Get taps by freq sampling
eqTaps = fir2(numTaps, linspace(0, 1, nfft/2+1), HinvWeighted.');
if MINPHASE
    [~, eqTaps] = rceps(eqTaps);
end
eqTaps = eqTaps.';
end
