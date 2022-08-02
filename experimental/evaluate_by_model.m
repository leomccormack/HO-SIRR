


%% Read
ref_order = 8;  % highest available

[irR, fs] = audioread(['MagLS_o', num2str(ref_order), '.wav']);

[irA1, fs] = audioread(['HOSIRR_o', num2str(3), '_sirr_vls.wav']);
[irA2, fs] = audioread(['HOSIRR_o', num2str(3), '_sirr_bind.wav']);
[irA3, fs] = audioread(['HOSIRR_o', num2str(5), '_sirr_vls.wav']);
[irA4, fs] = audioread(['HOSIRR_o', num2str(5), '_sirr_bind.wav']);
numSmps = size(irR,1);

%% create model input
s = randn(fs, 1);
s = cat(1, s, zeros(numSmps,1));
tsA1 = fftfilt(irA1, s);
tsA2 = fftfilt(irA2, s);
tsA3 = fftfilt(irA3, s);
tsA4 = fftfilt(irA4, s);

tsR = fftfilt(irR, s);

% combine stimuli into one matrix
ts = cat(3,tsA1,tsA2,tsA3,tsA4);
rs = cat(3,tsR,tsR,tsR,tsR);
tsP = permute(ts,[1 3 2]);
rsP = permute(rs,[1 3 2]);


%% CLL
[tsRbands, fc] = auditoryfilterbank(tsR, fs, 'fhigh', 20000);
[tsA1bands, fc] = auditoryfilterbank(tsA1, fs, 'fhigh', 20000);
[tsA2bands, fc] = auditoryfilterbank(tsA2, fs, 'fhigh', 20000);
[tsA3bands, fc] = auditoryfilterbank(tsA3, fs, 'fhigh', 20000);
[tsA4bands, fc] = auditoryfilterbank(tsA4, fs, 'fhigh', 20000);


CLL_R = calcCLL(tsRbands);
CLL_A1 = calcCLL(tsA1bands);
CLL_A2 = calcCLL(tsA2bands);
CLL_A3 = calcCLL(tsA3bands);
CLL_A4 = calcCLL(tsA4bands);

figure; hold on
plot(fc, CLL_A1 - CLL_R, 'LineWidth', 2);
plot(fc, CLL_A2 - CLL_R, 'LineWidth', 2);
plot(fc, CLL_A3 - CLL_R, '--', 'LineWidth', 2);
plot(fc, CLL_A4 - CLL_R, '--', 'LineWidth', 2);
set(gca(),'xscale','log')
title("CLL Difference")
grid on
legend(["VLS3", "DBIN3", "VLS5", "DBIN5"])
xlabel("Frequency in Hz")
ylabel("phons")
%figure; plot(db(CLL_A1 ./ CLL_R)); hold on;  plot(db(CLL_A3 ./ CLL_R))


%% Parameters
%amt_start

domFlag = 0; % specify that inputs are time-domain signals
freqRange = [100 15000]; % calculate spectral difference between 20Hz and 20kHz
nfft = length(rs(:,1,1)); % fft window size same as signal length
fpars.fs = fs;fpars.nfft = nfft;fpars.minFreq = freqRange(1); fpars.maxFreq = freqRange(2);
datasetNormalisation = 0; % blank vector for using iterative dataset normalisation. if an int, then that fixes the dataset normalisation in dB. Thus for no normalisation, set to 0.

% Calculate perceptual spectral difference
[specDiff,PSpecDiff] = mckenzie2022(tsP,rsP,domFlag,fpars,datasetNormalisation, 1, 1);
PSpecDiff = squeeze(PSpecDiff);

% get single values of spectral difference for all stimuli
%PavgSpecDiffS = mean(PSpecDiff,2);
disp('===')
disp(PSpecDiff)


%%
function CLL = calcCLL(x_bands)
CLL = log2(squeeze( sum(sqrt(sqrt(sum(x_bands.^2,1) / size(x_bands,1))), 3) )) *10+40;
end
