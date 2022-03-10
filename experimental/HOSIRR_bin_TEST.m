%%%%% 
clear all, close all, dbstop if error %#ok

addpath '..'
addpath '../_Simulated_Rooms_' '../_Stimuli_'

demo_order = 5;  
[input_stimulus, fs_in] = audioread('music__KickDrumClicky.wav');


%% FIRST-ORDER SIRR BIN CHECKS
pars.order = 1;  
pars.fs = 48e3; 
[~,dirs_rad] = getTdesign(2*(pars.order+1));
pars.ls_dirs_deg = dirs_rad*180/pi;
pars.multires_winsize = 128;  
pars.multires_xovers = [];   
pars.RENDER_DIFFUSE = 1;
pars.BROADBAND_FIRST_PEAK = 0; %%%DISABLED     
pars.nBroadbandPeaks = 1;    
pars.decorrelationType = 'noise'; 
pars.BROADBAND_DIFFUSENESS = 1;
pars.maxDiffFreq_Hz = 3000;  
pars.alpha_diff = 0.5;
pars.chOrdering = 'ACN';
pars.normScheme = 'N3D';
pars.hrtf_sofa_path = '/Users/mccorml1/Documents/HRIRs_SOFA/kemarhead_aalto2016.sofa';
pars.hrtf_sofa_path = '/home/chris/data/HRTFs/Kemar_Aalto_2016/kemarhead_aalto2016.sofa';

% --- Single plane-wave input --- 
src_dir = [-45 -45];  % try adding more than 1 plane-wave, to see the first-order analysis break
shir = randn(pars.fs, size(src_dir,1)) * (sqrt(4*pi).*getRSH(pars.order, src_dir)).';
[~,~,~,~,analysis] = HOSIRR_bin(shir, pars);

% In this case, we would expect:
% - [azimuth elevation] should correspond to 'src_dir'  
% - diffuseness should be 0
% - non-diffuse energy should be similar to input sound-field energy, and 
%   diffuse energy ~0
%figure, subplot(4,1,1), imagesc(analysis.azim{1}*180/pi), colorbar, axis xy, caxis([-180 180]), title('azimuth (degrees)')
%subplot(4,1,2), imagesc(analysis.elev{1}*180/pi), colorbar, axis xy, caxis([-90 90]), title('elevation (degrees)')
%subplot(4,1,3), imagesc(10*log10(analysis.energy{1})), colorbar, axis xy, title('energy (dB)')
%subplot(4,1,4), imagesc(analysis.diff{1}), colorbar, axis xy, caxis([0 1]), title('diffuseness')
figure, plot(10*log10(analysis.sf_energy{1})), hold on 
plot(10*log10(analysis.ndiff_energy{1})), hold on
plot(10*log10(analysis.diff_energy{1})) 
title('energy (dB)'), grid on, ylim([-40 20])
legend('sound-field', 'non-diffuse', 'diffuse')
%figure, plot(analysis.diff{1}(1,:)), title('diffuseness'), grid on, ylim([0 1])

% --- Diffuse input --- 
[~, diff_dirs] = getTdesign(21); % approximate diffuse-field with 240 incoherent noise sources
shir = randn(pars.fs, size(diff_dirs,1)) * (sqrt(4*pi).*getRSH(pars.order, diff_dirs*180/pi)).';
[~,~,~,~,analysis] = HOSIRR_bin(shir, pars);

% In this case, we would expect:
% - [azimuth elevation] should be random 
% - diffuseness should be close to 1
% - diffuse energy should be similar to the input sound-field energy, and
%   non-diffuse energy much lower than diffuse energy
%figure, subplot(4,1,1), imagesc(analysis.azim{1}*180/pi), colorbar, axis xy, caxis([-180 180]), title('azimuth (degrees)')
%subplot(4,1,2), imagesc(analysis.elev{1}*180/pi), colorbar, axis xy, caxis([-90 90]), title('elevation (degrees)')
%subplot(4,1,3), imagesc(10*log10(analysis.energy{1})), colorbar, axis xy, title('energy (dB)')
%subplot(4,1,4), imagesc(analysis.diff{1}), colorbar, axis xy, caxis([0 1]), title('diffuseness')
figure, plot(10*log10(analysis.sf_energy{1})), hold on 
plot(10*log10(analysis.ndiff_energy{1})), hold on
plot(10*log10(analysis.diff_energy{1})) 
title('energy (dB)'), grid on, ylim([-40 20])
legend('sound-field', 'non-diffuse', 'diffuse')
%figure, plot(analysis.diff{1}(1,:)), title('diffuseness'), grid on, ylim([0 1])

clear pars

%% EITHER: CREATE DEMO SPHERICAL HARMONIC (SH; AMBISONIC/B-FORMAT) RIR
% For demonstration purposes, this is a simulated Auditorium using LoRA [1]
% configured for DTU's AVIL anechoic listening room.
[ref_ls_rir, fs] = audioread('Auditorium_64ch_DTU_AVIL.wav');  
load('DTU_ls_dirs_deg.mat')

% Convolve the input stimulus with this reference loudspeaker array RIR.
% This will serve as the reference rendering.
assert(fs == fs_in)
ir_length = size(ref_ls_rir, 1);
nLS = size(ref_ls_rir, 2);
output = matrixConvolver(input_stimulus, reshape(ref_ls_rir, [ir_length, 1, nLS]), size(input_stimulus,1));
output = 0.99.*output./max(abs(output(:)));
audiowrite('BFormatDemo_reference.wav', output, fs, 'BitsPerSample', 24);

% encode into the spherical harmonic domain
sh_rir = ref_ls_rir * sqrt(4*pi/nLS) * getRSH(demo_order, ls_dirs_deg)'; 
sh_rir = 0.99.*sh_rir./max(abs(sh_rir(:)));
clear ls_dirs_deg;


%% OR: LOAD  ROOM IMPULSE RESPONSE FROM FILE 
%[sh_rir, fs] = audioread('PATH_TO_SH_RIR.wav');
disp(' * Replace "sh_rir" with your own Ambisonic/B-Formt room impulse response measurement')


%% DEFINE HOSIRR USER PARAMETERS
% Specify input signal channel ordering and normalisation conventions
% Note: 
%     Ambix -> ACN/SN3D
%     FuMa  -> WXYZ/SN3D (first-order only. Note this is also without the 
%     1/sqrt(2) scaling on the omni component)
pars.chOrdering = 'ACN'; % 'ACN', or 'WXYZ' (deprecated/first-order only)
pars.normScheme = 'N3D'; % 'N3D', or 'SN3D'
pars.fs = fs;  
% Specify windowing size, in samples, (note HOSIRR employs 50% overlap)
pars.multires_winsize = 128;
pars.multires_xovers = [ ];   
% or if you want to use the multi-resolution STFT option, e.g.:
%    512 samples hop size up to 500 Hz, then:
%    128 samples hop size   "   2 kHz,    "
%    64  samples hop size above 2 kHz
% then set the following:
%pars.multires_winsize = [512, 128, 64]; 
%pars.multires_xovers = [500, 2e3];   
pars.RENDER_DIFFUSE = 1;
pars.decorrelationType = 'noise';
% This option isolates the first peak in the response and renders it based
% on a broad-band DoA estimate. 
pars.BROADBAND_FIRST_PEAK = 1;  
% This option allows the diffuseness parameter to be computed broad-band
% (up to "maxDiffFreq_Hz") and replicated to all bin, rather than
% computing the diffuseness for each bin independently
pars.BROADBAND_DIFFUSENESS = 1;
pars.maxDiffFreq_Hz = 3000;  
% diffuseness parameter temporal averaging coefficient (one-pole filter)
pars.alpha_diff = 0.5; 

% 
% %% EITHER: DEFINE THE SAME LOUDSPEAKER DIRECTIONS AS REFERENCE
load('DTU_ls_dirs_deg.mat') 
[~, dirs_rad] = getTdesign(5);
pars.ls_dirs_deg = dirs_rad*180/pi;
% 
% 
% %% OR: DEFINE LOUDSPEAKER DIRECTIONS FOR YOUR SETUP
% %pars.ls_dirs_deg = [];
% disp(' * Replace "ls_dirs_deg" with your own loudspeaker array angles')
% % Define your loudspeaker directions in the [azim elev] convention, in 
% % DEGREES; nLoudspeakers x 2)
% % pars.ls_dirs_deg = [ 0 0; 45 0; -45 0; 90 0; -90 0; 135 0; -135 0; 45 35; -45 35; 90 35; -90 35; 0 90;]; 
% 
% 
% %% RENDER SH RIR TO TARGET LOUDSPEAKER LAYOUT
% % Render
tic, sirr_ls_rir = HOSIRR(sh_rir, pars); toc
audiowrite(['HOSIRR_ls_o' num2str(demo_order) '.wav'], 0.9.*sirr_ls_rir, fs);
% 
% % Convolve the input stimulus with the rendered loudspeaker array RIR
% ir_length = size(sirr_ls_rir, 1);
% nLS = size(sirr_ls_rir, 2);  
% output = matrixConvolver(input_stimulus, reshape(sirr_ls_rir, [ir_length, 1, nLS]), size(input_stimulus,1));
% output = 0.99.*output./max(abs(output(:))); % normalise
% audiowrite(['BFormatDemo_HOSIRR_o' num2str(demo_order) '.wav'], output, fs, 'BitsPerSample', 24);

%% WIP: BINAURAL
disp(' * BINAURAL')
% load HRIRs
%pars.hrtf_sofa_path = '/Users/holdc1/Documents/data/HRTFs/Kemar_Aalto_2016/kemarhead_aalto2016.sofa';
%pars.hrtf_sofa_path = '/Users/mccorml1/Documents/HRIRs_SOFA/kemarhead_aalto2016.sofa';
pars.hrtf_sofa_path = '/home/chris/data/HRTFs/Kemar_Aalto_2016/kemarhead_aalto2016.sofa';

assert(isfile(pars.hrtf_sofa_path))
pars.BROADBAND_FIRST_PEAK = 0;
[sirr_bin, ~, ~, pars] = HOSIRR_bin(sh_rir, pars);
audiowrite(['HOSIRR_o' num2str(demo_order) '_bin.wav'], 0.9.*sirr_bin, fs);

LISTEN = true
if LISTEN
hrir_0 = pars.hrirs(:, :, 6);
mono_bir = fftfilt(hrir_0, sqrt(4*pi)*sh_rir(:, 1));
sound(mono_bir, fs)
pause(2)
sound(sirr_bin, fs)
pause(2)
end

figure
subplot(2,1,1)
plot(mono_bir)
subplot(2,1, 2)
plot(sirr_bin)
rms(mono_bir(0.1*fs:end-100, :))
rms(sirr_bin(0.1*fs:end-100, :))


%% Compare
[x_ls, y_ls, z_ls] = sph2cart( pars.ls_dirs_deg(:,1)*pi/180,...
                               pars.ls_dirs_deg(:,2)*pi/180, 1);
[x_hrfts, y_hrtfs, z_hrtfs] = sph2cart(pars.hrtf_dirs_deg(:, 1)*pi/180,...
                                       pars.hrtf_dirs_deg(:, 2)*pi/180, 1);
ls_proj = [x_hrfts, y_hrtfs, z_hrtfs] * [x_ls, y_ls, z_ls].';
[d_min, d_min_k] = max(ls_proj);
ls_sigs_bin_l = sum(fftfilt(squeeze(pars.hrirs(:, 1, d_min_k)), sirr_ls_rir), 2);
ls_sigs_bin_r = sum(fftfilt(squeeze(pars.hrirs(:, 2, d_min_k)), sirr_ls_rir), 2);
disp("LS HOSIRR")
sound([ls_sigs_bin_l, ls_sigs_bin_r], fs)
pause(2)
disp("Binaural HOSIRR")
sound(sirr_bin, fs)
pause(2)

% Prepare HRTFs
hrtfs = fft(pars.hrirs, [], 1);
%hrtfs_syn = fft(hrir, fftsize, 1);
if mod(size(pars.hrirs, 1), 2)
    error('not implemented')
else  % even
    hrtfs = hrtfs(1:size(pars.hrirs, 1)/2+1, :, :);  % pos half
end
D_bin = getAmbisonic2BinauralFilters_magls_zotter(permute(hrtfs, [2 3 1]),...
    pars.hrtf_dirs_deg, pars.order, [], pars.fs, pars.hrirs_weights);
D_bin = real(ifft(cat(3, D_bin, conj(D_bin(:,:, end-1:-1:2))), [], 1));
sh_rir_bin_l = sum(fftfilt(squeeze(D_bin(1, :,:)).', sh_rir), 2);
sh_rir_bin_r = sum(fftfilt(squeeze(D_bin(2, :,:)).', sh_rir), 2);
disp("MagLS")
soundsc([sh_rir_bin_l, sh_rir_bin_r], fs)
pause(2)
%out2 = matrixConvolver(sh_rir, permute(D_bin,[3, 2, 1]), 2048);
%soundsc(out2, fs)

