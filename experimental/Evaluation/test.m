clear all, close all %#ok
% first run: "ims_test_IRs/enerate_test_irs.m" to generate test RIRs

addpath('../')           % for HOSIRR_bin.m 
addpath('../../')        % for matrixconv.m 
addpath('ims_test_IRs/') % for the reference RIRs
addpath('HighTdesigns/') % for the evaluation grid

% Config
%ref_name = 'ref_o20_medium_room_leo';
ref_name = 'ref_o20_medium_room';

% Load reference RIR
[ref_rir, fs] = audioread([ref_name '.wav']);
ref_order = 5;

% Pass order truncated reference through HOSIRR_bin
trunc_order = 1;
pars.chOrdering = 'ACN'; 
pars.normScheme = 'N3D'; 
pars.fs = fs;   
pars.BROADBAND_FIRST_PEAK = 0;
pars.RENDER_DIFFUSE = 1;
pars.decorrelationType = 'noise';
pars.BROADBAND_DIFFUSENESS = 1;
pars.maxDiffFreq_Hz = 8000;   
pars.alpha_diff = 0.5;
pars.multires_winsize = 128;  
pars.multires_xovers = [ ];      
pars.hrtf_sofa_path = '/Users/mccorml1/Documents/HRIRs_SOFA/kemarhead_aalto2016.sofa';
tic, [sirr_bin_rir,~,~,~,analysis] = HOSIRR_bin(ref_rir(:,1:(trunc_order+1).^2), pars); toc
sirr_bin_rir = sirr_bin_rir./max(abs(sirr_bin_rir(:)));
audiowrite('sirr_bin_rir.wav', sirr_bin_rir, fs, 'BitsPerSample', 24);
 
% MagLS for comparison
hrirs = ncread(pars.hrtf_sofa_path,'Data.IR');
hrir_dirs_deg = ncread(pars.hrtf_sofa_path,'SourcePosition'); 
hrtf_dirs_deg = hrir_dirs_deg(1:2,:).'; 
hrtf_dirs_rad = hrtf_dirs_deg*pi/180;
hrtf_dirs_rad2 = [hrtf_dirs_rad(:,1) pi/2-hrtf_dirs_rad(:,2)];
nHRTF = size(hrtf_dirs_deg,1);
lHRTF = size(hrirs,1);
%w = (getVoronoiWeights(hrtf_dirs_rad));
w = ones(size(hrtf_dirs_rad,1),1)./size(hrtf_dirs_rad,1);
hrtfs = fft(hrirs,[],1);  
hrtfs = hrtfs(1:lHRTF/2+1,:,:);
hrtfs = permute(hrtfs,[2 3 1]);
f = (0:lHRTF/2)'*fs/lHRTF;
[D_magls, h_magls] = getAmbisonic2BinauralFilters_magls_zotter(hrtfs, hrtf_dirs_deg, trunc_order, [], fs, w);
[D_ref, h_ref] = getAmbisonic2BinauralFilters_magls_zotter(hrtfs, hrtf_dirs_deg, ref_order, [], fs, w);

% MagLS at same order as HOSIRR_bin
magls_bin_rir = matrixConvolver(ref_rir(:,1:(trunc_order+1).^2), h_magls, size(h_magls,1));
magls_bin_rir = magls_bin_rir./max(abs(magls_bin_rir(:)));
audiowrite('magls_bin_rir.wav', magls_bin_rir, fs, 'BitsPerSample', 24);

% "Reference" being higher order MagLS, the current HRTF set is not dense enough to
% go really high here
ref_bin_rir = matrixConvolver(ref_rir(:,1:(ref_order+1).^2), h_ref, size(h_ref,1));
ref_bin_rir = ref_bin_rir./max(abs(ref_bin_rir(:)));
audiowrite('ref_bin_rir.wav', ref_bin_rir, fs, 'BitsPerSample', 24);



