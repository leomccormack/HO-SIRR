function [pars] = prepareHRIRs(pars, sofapath)
%UNTITLED Summary of this function goes here
%   Detailed explanation goes here
assert(isfile(sofapath), "SOFA file not found!")
pars.hrtf_sofa_path = sofapath;
hrirs = ncread(pars.hrtf_sofa_path,'Data.IR');
hrir_dirs_deg = ncread(pars.hrtf_sofa_path,'SourcePosition');
pars.numHrirs = size(hrirs, 3);
pars.lenHrirs = size(hrirs, 1);


hrir_dirs_deg = hrir_dirs_deg(1:2,:).'; %         % nHRTF x 2     (deg)
pars.hrir_dirs_deg = hrir_dirs_deg;
%pars.hrirs_weights = getVoronoiWeights(pars.hrtf_dirs_deg);
[pars.hrirs_weights,N_support] = findGridWeights(deg2rad(hrir_dirs_deg(:,1)),...
                                     pi/2-deg2rad(hrir_dirs_deg(:,2)));

itd = computeITDfromXCorr(hrirs, pars.fs);
pars.hrtf_itd = itd;        % nHRTF x 1
eqTaps = hrirsDiffuseFieldEQ(hrirs, 1, pars.hrirs_weights);
tapswin = cos(linspace(0, pi/2, 32));
eqTaps = eqTaps(1:256);
eqTaps(end-31:end) = eqTaps(end-31:end) .* tapswin.';
hrirs_filt(:,1,:) = fftfilt(eqTaps, ...
    cat(1,squeeze(hrirs(:,1,:)),zeros(length(eqTaps),pars.numHrirs)));
hrirs_filt(:,2,:) = fftfilt(eqTaps, ...
    cat(1,squeeze(hrirs(:,2,:)),zeros(length(eqTaps),pars.numHrirs)));

% write back
pars.hrirs = hrirs_filt;
end