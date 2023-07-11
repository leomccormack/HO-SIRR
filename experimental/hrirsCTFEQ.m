function [eqTaps] = hrirsCTFEQ(hrirs, fs, MINPHASE, gridWeights)
%hrirsCTFEQ Calculate diffuse field (common transfer function) EQ.
%   INPUTS
%       hrirs : [len x 2 x grid]
%       MINPHASE : bool (optional, default true)
%       gridWeights : [grid, 1] (optional, default [] assumes regular grid)
%   OUTPUT
%       eqTaps : Filter taps [len x 1]
%       
% Chris Hold 2021,2022,2023

if nargin <4
    gridWeights = [];
end
if nargin < 3
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
nfft = max(2048, 16*numTaps);  % interpolate
H = fft(hrirs, nfft, 1);
Hs = H(1:nfft/2+1, :, :);


% weighted RMS
Havg = sqrt(sum(reshape(gridWeights, 1, 1, []) .* abs(Hs).^2, 3) / (4*pi));

% Smoothing
HavgSmooth = Havg;
for i = 1:2
for bin = 2:nfft/2+1
    if ~mod(bin, 2)
        avgidx = (bin - bin/2 : min(3/2 * bin, nfft/2+1));
    end
    win = hann(length(avgidx));
    % weighted average
    HavgSmooth(bin, i) = sum(win .* Havg(avgidx, i)) / sum(win);
end
end

% Avg (left, right)
HavgSmooth = mean(Havg, 2);


% frequency mask for lo and hi, avoid inversion there
freqWeight = ones(nfft/2+1, 1);
fv = linspace(0, fs/2, nfft/2+1);
fLo = 100;
[~, midx] = min(abs(fv - fLo));
wLo = hann(2*midx+1);
freqWeight(1:midx + 1) = wLo(1:midx + 1);
fHi = 15000;
[~, midx] = min(abs(fv - fHi));
midx = length(fv) - midx;
wHi = hann(2*midx+1);
freqWeight(end-midx:end) = wHi(end-midx:end);

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
