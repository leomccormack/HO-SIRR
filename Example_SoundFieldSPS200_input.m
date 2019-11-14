%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% SoundField SPS200 input example
% -------------------------------
% An example of how to render input SoundField SPS200 room impulse 
% responses (RIRs) for a specific loudspeaker arrangement using HO-SIRR.
% Note: you may alter the array radius and/or the mic channel order, in
% order to cater for different A-Format arrays with the same script.
% 
% DEPENDENCES
%   Spherical-Harmonic-Transform Matlab library
%       https://github.com/polarch/Spherical-Harmonic-Transform
%   Higher-Order-Ambisonics Matlab library
%       https://github.com/polarch/Higher-Order-Ambisonics
%   Array-Response-Simulator
%       https://github.com/polarch/Array-Response-Simulator
%   Vector-Base-Amplitude-Panning
%       https://github.com/polarch/Vector-Base-Amplitude-Panning
%
% REFERENCES
%   [1] Favrot, S. and Buchholz, J.M., 2010. 
%       LoRA: A loudspeaker-based room auralization system. Acta Acustica 
%       united with Acustica, 96(2), pp.364-375.
%   [2] Moreau, S., Daniel, J., Bertet, S., 2006, 
%       3D sound field recording with higher order ambisonics-objective 
%       measurements and validation of spherical microphone. In Audio 
%       Engineering Society Convention 120.
%   [3] Politis, A., Gamper, H. 2017. 
%       "Comparing Modelled And Measurement-Based Spherical Harmonic 
%       Encoding Filters For Spherical Microphone Arrays. In IEEE 
%       Workshop on Applications of Signal Processing to Audio and 
%       Acoustics (WASPAA).
%   [4] sparta_binauraliser VST plug-in: 
%       https://github.com/leomccormack/SPARTA
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
%   Leo McCormack, 28/10/2019
%   leo.mccormack@aalto.fi 
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
clear all, close all, dbstop if error %#ok
addpath '_Simulated_Rooms_' '_Stimuli_'

demo_order = 1; 
nSH = int32((demo_order+1).^2);
[input_stimulus, fs_in] = audioread('music__KickDrumClicky.wav');

% SoundField SPS200 angles, in [azimuth, elevation]:
mic_dirs_deg = ... 
    [45 35.264; -45 -35.264; 135 -35.264; -135 35.264;];
mic_dirs_rad = mic_dirs_deg*pi/180;
R = 0.02; % the SoundField SPS200 radius 
arrayType = 'directional'; % the SoundField SPS200 is an open directional (cardioid) array


%% EITHER: CREATE DEMO SOUNDFIELD SPS200 ROOM IMPULSE RESPONSE (RIR)
% This is a simulated Auditorium room using ODEON and LoRA [1], which is
% configured for DTU's AVIL anechoic listening room. 
[ref_ls_rir, fs] = audioread('Auditorium_64ch_DTU_AVIL.wav');  
load('DTU_ls_dirs_deg.mat') 
ls_dirs_rad = ls_dirs_deg*pi/180; % loudspeaker directions for the AVIL setup

% Convolve the input stimulus with this reference loudspeaker array RIR.
% This will serve as the reference rendering.
assert(fs == fs_in)
ir_length = size(ref_ls_rir, 1);
nLS = size(ref_ls_rir, 2);
output = matrixConvolver(input_stimulus, reshape(ref_ls_rir, [ir_length, 1, nLS]), size(input_stimulus,1));
output = 0.99.*output./max(abs(output(:)));
audiowrite('SoundFieldSPS200Demo_reference.wav', output, fs, 'BitsPerSample', 24);

% We then simulate an SPS200 array, deriving the theoretical array transfer 
% functions for each loudspeaker direction in the DTU set-up
N_order = 15; % order of approximation 
Lfilt = 1024; % filter length
[h_mic, H_mic] = simulateSphArray(Lfilt, mic_dirs_rad, ls_dirs_rad, arrayType, R, N_order, fs, 0.5);

% By convolving each loudspeaker RIR with the corresponding SPS200
% transfer functions, we essentially simulate an SPS200 RIR 
% "measurement" of the reference room
mic_rir = matrixConvolver(ref_ls_rir, permute(h_mic, [1 3 2]), Lfilt/2);
mic_rir = 0.99.*mic_rir./max(abs(mic_rir(:)));

clear ls_dirs_deg ls_dirs_rad


%% OR: LOAD SOUNDFIELD SPS200 RIR FROM FILE 
%[mic_rir, fs] = audioread('PATH_TO_SPS200_RIR.wav');
disp(' * Replace "mic_rir" with your own SoundField SPS200 room impulse response measurement')


%% SOUNDFIELD SPS200 RIR TO SPHERICAL HARMONIC RIR CONVERSION
Lfilt = 1024; % filter length
maxG_dB = 10; % maximum allowed gain amplification of the encoding filters
nMics = size(mic_dirs_rad,1);

% modal responses
c = 343;
f = (0:Lfilt/2)'*fs/Lfilt;
kR = 2*pi*f*R/c;
bN = sphModalCoeffs(demo_order, kR, arrayType, 0.5)/(4*pi);

% encoding matrix with regularization 
alpha = sqrt(nMics)*10^(maxG_dB/20);
beta = sqrt((1-sqrt(1-1/alpha^2))/(1+sqrt(1-1/alpha^2)))+eps; % Moreau & Daniel
% regularized single channel equalization filters per order
H_filt = conj(bN)./(abs(bN).^2 + beta^2*ones(Lfilt/2+1,demo_order+1));

% time domain filters
h_filt = H_filt;
h_filt(end,:) = real(h_filt(end,:));
h_filt = [h_filt; conj(h_filt(end-1:-1:2,:))];
h_filt = real(ifft(h_filt));
h_filt = fftshift(h_filt, 1);    

% plot encoding filters
figure
subplot(2,1,1)
bN_inv_dB = 20*log10(abs(H_filt(2:end,:)));
semilogx(f(2:end,1).',bN_inv_dB)
grid on
set(gca,'xlim',[20 20000],'ylim',[min(min(bN_inv_dB))-10 max(max(bN_inv_dB))+10])
xlabel('\omega (Hz)'), ylabel('|H(\omega)| (dB)') 
subplot(2,1,2)
bN_dB = 20*log10(abs(1./(bN(2:end,:).'+eps)));
semilogx(f(2:end,1).',bN_dB)
grid on 
set(gca,'xlim',[20 20000],'ylim',[min(min(bN_dB))-10 max(max(bN_dB))+10])
xlabel('\omega (Hz)'), ylabel('|H(\omega)| (dB)')
title(' no reg') 

% get real SH matrix for microphone directions and replicate encoding
% filters for each order
aziElev2aziPolar = @(dirs) [dirs(:,1) pi/2-dirs(:,2)]; % macro to convert from azimuth-inclination to azimuth-elevation
Y_mics = sqrt(4*pi) * getSH(demo_order, aziElev2aziPolar(mic_dirs_rad), 'real'); 
h_filt_r = replicatePerOrder(h_filt/sqrt(4*pi), 2);

% apply spherial harmonic transform
sh_rir_tmp = mic_rir * pinv(Y_mics).';
sh_rir = zeros(size(sh_rir_tmp,1)+Lfilt-1, nSH);
for i=1:nSH
    sh_rir(:,i) = fftconv(sh_rir_tmp(:,i), h_filt_r(:,i));  
end
sh_rir = 0.99.*sh_rir./max(abs(sh_rir(:)));
disp(' * Or, replace "sh_rir" with an encoding of your own SoundField SPS200 room impulse response measurement')

 
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
% (up to "maxDiffFreq_Hz") and replicated to all bins, rather than
% computing the diffuseness for each bin independently
pars.BROADBAND_DIFFUSENESS = 1;
pars.maxDiffFreq_Hz = 3000;  
% diffuseness parameter temporal averaging coefficient (one-pole filter)
pars.alpha_diff = 0.5;


%% EITHER: DEFINE THE SAME LOUDSPEAKER DIRECTIONS AS REFERENCE
load('DTU_ls_dirs_deg.mat') 
pars.ls_dirs_deg = ls_dirs_deg;


%% OR: DEFINE LOUDSPEAKER DIRECTIONS FOR YOUR SETUP
%pars.ls_dirs_deg = [];
disp(' * Replace "ls_dirs_deg" with your own loudspeaker array angles')
% Define your loudspeaker directions in the [azim elev] convention, in 
% DEGREES; nLoudspeakers x 2)


%% RENDER SH RIR TO TARGET LOUDSPEAKER LAYOUT
tic, sirr_ls_rir = HOSIRR(sh_rir, pars); toc 
%audiowrite(['HOSIRR_o' num2str(pars.maxOrder) '.wav'], 0.9.*sirr_ls_rir, fs);  

% Convolve the input stimulus with the rendered loudspeaker array RIR
ir_length = size(sirr_ls_rir, 1);
nLS = size(sirr_ls_rir, 2);  
output = matrixConvolver(input_stimulus, reshape(sirr_ls_rir, [ir_length, 1, nLS]), size(input_stimulus,1));
output = 0.99.*output./max(abs(output(:))); % normalise
audiowrite(['SoundFieldSPS200Demo_HOSIRR_o' num2str(demo_order) '.wav'], output, fs, 'BitsPerSample', 24);

% The reference and HOSIRR renderings may then be compared over headphones using
% e.g. the free sparta_binauraliser VST plug-in [4] (which has a DTU setup preset).  
% Or, of course, you may also compare the two renders using the actual DTU 
% loudspeaker setup :-)
