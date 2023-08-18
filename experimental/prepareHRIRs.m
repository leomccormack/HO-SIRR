function [pars] = prepareHRIRs(pars, sofapath)
%UNTITLED Summary of this function goes here
%   Detailed explanation goes here
if isfield(pars, 'hrirs'); warning("HRIRs already prepared, overwriting."); end
assert(isfile(sofapath), "SOFA file not found!")
pars.hrtf_sofa_path = sofapath;
hrirs = ncread(pars.hrtf_sofa_path,'Data.IR');
hrir_dirs_deg = ncread(pars.hrtf_sofa_path,'SourcePosition');
fs_sofa = ncread(pars.hrtf_sofa_path,'Data.SamplingRate');
assert(fs_sofa == pars.fs)
numHrirs = size(hrirs, 3);

hrir_dirs_deg = hrir_dirs_deg(1:2,:).'; %         % nHRTF x 2     (deg)
pars.hrir_dirs_deg = hrir_dirs_deg;
%pars.hrirs_weights = getVoronoiWeights(pars.hrtf_dirs_deg);
[pars.hrirs_weights,N_support] = findGridWeights(deg2rad(hrir_dirs_deg(:,1)),...
                                     pi/2-deg2rad(hrir_dirs_deg(:,2)));

itd = computeITDfromXCorr(hrirs, pars.fs);
pars.hrtf_itd = itd;        % nHRTF x 1
eqTaps = hrirsCTFEQ(hrirs, pars.fs, 0, pars.hrirs_weights);
eqTaps = eqTaps(1:size(hrirs,1));
nTaperTaps = round(size(hrirs,1) / 8);
tapswin = cos(linspace(0, pi/2, nTaperTaps));
eqTaps(end-nTaperTaps+1:end) = eqTaps(end-nTaperTaps+1:end) .* tapswin.';
hrirs_filt(:,1,:) = fftfilt(eqTaps, ...
    cat(1,squeeze(hrirs(:,1,:)),zeros(length(eqTaps), numHrirs)));
hrirs_filt(:,2,:) = fftfilt(eqTaps, ...
    cat(1,squeeze(hrirs(:,2,:)),zeros(length(eqTaps), numHrirs)));

% write back
pars.hrirs = hrirs_filt;
pars.numHrirs = size(pars.hrirs, 3);
pars.lenHrirs = size(pars.hrirs, 1);
pars.CTFeqTaps = eqTaps;
end