%%%%% 
clear all, close all %, dbstop if error %#ok

addpath '..'
addpath '../_Simulated_Rooms_' '../_Stimuli_'

demo_order = 3;  
[input_stimulus, fs_in] = audioread('music__KickDrumClicky.wav');
%sofapath = '~/data/HRTFs/Kemar_Aalto_2016/KemarAuralID.sofa';
sofapath = '~/data/HRTFs/THK_KU100/HRIR_L2354.sofa';


% %% FIRST-ORDER SIRR BIN CHECKS
% pars.order = 1;  
% pars.fs = 48e3; 
% [~,dirs_rad] = getTdesign(2*(pars.order+1));
% pars.ls_dirs_deg = dirs_rad*180/pi;
% pars.multires_winsize = 128;  
% pars.multires_xovers = [];   
% pars.RENDER_DIFFUSE = 1;
% pars.BROADBAND_FIRST_PEAK = 0; %%%DISABLED     
% pars.nBroadbandPeaks = 1;    
% pars.decorrelationType = 'phase';  %%%PHASE
% pars.BROADBAND_DIFFUSENESS = 0;
% pars.maxDiffFreq_Hz = 3000;  
% pars.alpha_diff = 0.5;
% pars.chOrdering = 'ACN';
% pars.normScheme = 'N3D';
% pars.hrtf_sofa_path = sofapath;
% 
% %% --- Single plane-wave input --- 
% src_dir = [-45 -45];  % try adding more than 1 plane-wave, to see the first-order analysis break
% shir = randn(pars.fs, size(src_dir,1)) * (sqrt(4*pi).*getRSH(pars.order, src_dir)).';
% [sirr,~,~,pars,analysis] = HOSIRR_bin(shir, pars);
% 
% % In this case, we would expect:
% % - [azimuth elevation] should correspond to 'src_dir'  
% % - diffuseness should be 0
% % - non-diffuse energy should be similar to input sound-field energy, and 
% %   diffuse energy ~0
% %figure, subplot(4,1,1), imagesc(analysis.azim{1}*180/pi), colorbar, axis xy, caxis([-180 180]), title('azimuth (degrees)')
% %subplot(4,1,2), imagesc(analysis.elev{1}*180/pi), colorbar, axis xy, caxis([-90 90]), title('elevation (degrees)')
% %subplot(4,1,3), imagesc(10*log10(analysis.energy{1})), colorbar, axis xy, title('energy (dB)')
% %subplot(4,1,4), imagesc(analysis.diff{1}), colorbar, axis xy, caxis([0 1]), title('diffuseness')
% figure, plot(10*log10(analysis.sf_energy{1})), hold on 
% plot(10*log10(analysis.ndiff_energy{1})), hold on
% plot(10*log10(analysis.diff_energy{1})) 
% title('PW: energy (dB)'), grid on
% legend('sound-field', 'non-diffuse', 'diffuse')
% %figure, plot(analysis.diff{1}(1,:)), title('diffuseness'), grid on, ylim([0 1])
% 
% %% --- Diffuse input --- 
% %[~, diff_dirs] = getTdesign(21); % approximate diffuse-field with 240 incoherent noise sources
% %shir = randn(pars.fs, size(diff_dirs,1)) * (sqrt(4*pi).*getRSH(pars.order, diff_dirs*180/pi)).';
% shir = randn(pars.fs, (pars.order+1)^2);
% [sirr,~,~,pars,analysis] = HOSIRR_bin(shir, pars);
% 
% % In this case, we would expect:
% % - [azimuth elevation] should be random 
% % - diffuseness should be close to 1
% % - diffuse energy should be similar to the input sound-field energy, and
% %   non-diffuse energy much lower than diffuse energy
% %figure, subplot(4,1,1), imagesc(analysis.azim{1}*180/pi), colorbar, axis xy, caxis([-180 180]), title('azimuth (degrees)')
% %subplot(4,1,2), imagesc(analysis.elev{1}*180/pi), colorbar, axis xy, caxis([-90 90]), title('elevation (degrees)')
% %subplot(4,1,3), imagesc(10*log10(analysis.energy{1})), colorbar, axis xy, title('energy (dB)')
% %subplot(4,1,4), imagesc(analysis.diff{1}), colorbar, axis xy, caxis([0 1]), title('diffuseness')
% figure, plot(10*log10(analysis.sf_energy{1})), hold on 
% plot(10*log10(analysis.ndiff_energy{1})), hold on
% plot(10*log10(analysis.diff_energy{1})) 
% title('Diffuse: energy (dB)'), grid on
% legend('sound-field', 'non-diffuse', 'diffuse')
% %figure, plot(analysis.diff{1}(1,:)), title('diffuseness'), grid on, ylim([0 1])

%%
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
% then set the following:
pars.RENDER_DIFFUSE = 1;
pars.decorrelationType = 'phase';
% This option isolates the first peak in the response and renders it based
% on a broad-band DoA estimate. 
pars.BROADBAND_FIRST_PEAK = 1;  
% This option allows the diffuseness parameter to be computed broad-band
% (up to "maxDiffFreq_Hz") and replicated to all bin, rather than
% computing the diffuseness for each bin independently
pars.BROADBAND_DIFFUSENESS = 0;
pars.maxDiffFreq_Hz = 3000;  
% diffuseness parameter temporal averaging coefficient (one-pole filter)
pars.alpha_diff = 0.5; 
% Sector pattern
pars.pattern = 'maxRE';

% 
% %% EITHER: DEFINE THE SAME LOUDSPEAKER DIRECTIONS AS REFERENCE
%load('DTU_ls_dirs_deg.mat') 
[~, dirs_rad] = getTdesign(2*demo_order+1);  % +1 ?
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
%% RENDER SH RIR TO TARGET LOUDSPEAKER LAYOUT
% % Render
tic, sirr_ls_rir = HOSIRR(sh_rir, pars); toc
%audiowrite(['HOSIRR_ls_o' num2str(demo_order) '.wav'], 0.9.*sirr_ls_rir, fs);
% 
% % Convolve the input stimulus with the rendered loudspeaker array RIR
% ir_length = size(sirr_ls_rir, 1);
% nLS = size(sirr_ls_rir, 2);  
% output = matrixConvolver(input_stimulus, reshape(sirr_ls_rir, [ir_length, 1, nLS]), size(input_stimulus,1));
% output = 0.99.*output./max(abs(output(:))); % normalise
% audiowrite(['BFormatDemo_HOSIRR_o' num2str(demo_order) '.wav'], output, fs, 'BitsPerSample', 24);

%% WIP: BINAURAL by Reencoding
disp(' * BINAURAL Reenc')
% load HRIRs
pars.hrtf_sofa_path = sofapath;

assert(isfile(pars.hrtf_sofa_path))
[sirr_binre, sir_ndiff, sir_diff, pars, analysis] = HOSIRR_bin_re(sh_rir, pars);
%audiowrite(['HOSIRR_o' num2str(demo_order) '_binre.wav'], 0.9.*sirr_binre, fs);

%hrir_0 = pars.hrirs(:, :, 6);
%omni_bir = fftfilt(hrir_0, sqrt(4*pi)*sh_rir(:, 1));
LISTEN = true
if LISTEN
%sound(omni_bir, fs)
%pause(2)
sound(sirr_binre, fs)
pause(2)
sound(sir_ndiff, fs)
pause(2)
sound(sir_diff, fs)
pause(2)
end


%% WIP: BINAURAL Direct
disp(' * BINAURAL Dir')
% load HRIRs
pars.hrtf_sofa_path = sofapath;

assert(isfile(pars.hrtf_sofa_path))
[sirr_bind, sir_ndiffd, sir_diffd, pars, analysis] = HOSIRR_bin_direct(sh_rir, pars);
%audiowrite(['HOSIRR_o' num2str(demo_order) '_bind.wav'], 0.9.*sirr_binre, fs);

%hrir_0 = pars.hrirs(:, :, 6);
%omni_bir = fftfilt(hrir_0, sqrt(4*pi)*sh_rir(:, 1));
LISTEN = true
if LISTEN
%sound(omni_bir, fs)
%pause(2)
sound(sirr_bind, fs)
pause(2)
sound(sir_ndiffd, fs)
pause(2)
sound(sir_diffd, fs)
pause(2)
end


%% Binauralize
hrtfs = fft(pars.hrirs, [], 1);
assert (mod(size(pars.hrirs, 1)+1, 2))
hrtfs = hrtfs(1:size(pars.hrirs, 1)/2+1, :, :);  % pos half

[D_bin, D_bin_filters] = getAmbisonic2BinauralFilters_magls_hfcont(permute(hrtfs, [2 3 1]),...
    pars.hrir_dirs_deg, pars.order, pars.fs, 1500, pars.hrirs_weights);
sh_rir_bin_l = sum(fftfilt(squeeze(D_bin_filters(1, :,:)).', sh_rir), 2);
sh_rir_bin_r = sum(fftfilt(squeeze(D_bin_filters(2, :,:)).', sh_rir), 2);
rir_magLS = [sh_rir_bin_l, sh_rir_bin_r];

[x_ls, y_ls, z_ls] = sph2cart( pars.ls_dirs_deg(:,1)*pi/180,...
                               pars.ls_dirs_deg(:,2)*pi/180, 1);
[x_hrfts, y_hrtfs, z_hrtfs] = sph2cart(pars.hrir_dirs_deg(:, 1)*pi/180,...
                                       pars.hrir_dirs_deg(:, 2)*pi/180, 1);
ls_proj = [x_hrfts, y_hrtfs, z_hrtfs] * [x_ls, y_ls, z_ls].';
[d_min, d_min_k] = max(ls_proj);
ls_sigs_bin_l = sum(fftfilt(squeeze(pars.hrirs(:, 1, d_min_k)), sirr_ls_rir), 2);
ls_sigs_bin_r = sum(fftfilt(squeeze(pars.hrirs(:, 2, d_min_k)), sirr_ls_rir), 2);
sirr_vls = [ls_sigs_bin_l, ls_sigs_bin_r];

%% Compare
if LISTEN
disp("MagLS")
sound(rir_magLS, fs)
pause(2)
disp("VLS HOSIRR")
sound(sirr_vls, fs)
pause(2)
disp("Binaural Reenc HOSIRR")
sound(sirr_binre, fs)
pause(2)
disp("Binaural Direct HOSIRR")
sound(sirr_bind, fs)
pause(2)
end
audiowrite(['MagLS_o', num2str(demo_order), '.wav'], rir_magLS, fs);
audiowrite(['HOSIRR_o', num2str(demo_order), '_sirr_vls.wav'], sirr_vls, fs);
audiowrite(['HOSIRR_o', num2str(demo_order), '_sirr_binre.wav'], sirr_binre, fs);
audiowrite(['HOSIRR_o', num2str(demo_order), '_sirr_bind.wav'], sirr_bind, fs);

%out2 = matrixConvolver(sh_rir, permute(D_bin_filters,[3, 2, 1]), 2048);
%soundsc(out2, fs)

% check tail rms, should be very close
%rms(omni_bir(0.5*fs:fs, :))
rms(rir_magLS(0.5*fs:fs, :))
rms(sirr_vls(0.5*fs:fs, :))
rms(sirr_binre(0.5*fs:fs, :))
rms(sirr_bind(0.5*fs:fs, :))

%%
figure
subplot(4,1,1)
plot(sirr_vls)
ylim([-1, 1])
grid on
title('VLS HOSIRR')
subplot(4,1, 2)
plot(sirr_binre)
ylim([-1, 1])
grid on
title('BIN SIRR')
subplot(4, 1, 3)
plot(sirr_bind)
ylim([-1, 1])
grid on
title('DBIN SIRR')
subplot(4,1,4)
plot(rir_magLS)
ylim([-1, 1])
grid on
title('MagLS')

%%
fidx = 5
numSecs = size(pars.sectorDirs, 1);
pdirs = cat(4,analysis.azim{1}, analysis.elev{1});
ps = (analysis.energy{1});
%pscale = 100;  % scale
pscale = 250*1/max(max(abs(ps(fidx, :,:))));
ps = pscale * ps;
ps(ps<10e-6) = 10e-6; 
pa = 1-(analysis.diff{1});
pa(pa<0.9) = 0.75 * pa(pa<0.9); 
figure
hold on
for idxs = 1:numSecs
s = scatter(rad2deg(squeeze(pdirs(fidx,:,idxs,1))), rad2deg(squeeze(pdirs(fidx,:,idxs,2))),...
            squeeze(ps(fidx,:,idxs)), 'filled');
s.AlphaData = squeeze(pa(fidx,:,idxs));
s.MarkerFaceAlpha = 'flat';
end
plot(rad2deg(pars.sectorDirs(:,1)), rad2deg(pars.sectorDirs(:,2)), 'k+', 'LineWidth',2)
alim([0.01, 1])

set(gca, 'XDir','reverse')
xticks([-180:90:180])
xlim([-180, 180])
yticks([-90:45:90])
ylim([-90, 90])
xlabel('Azimuth')
ylabel('Elevation')
grid on
daspect([1 1 1])
title("DOA (f=" + num2str(pars.centerfreqs_anl(fidx)) + "Hz)")

