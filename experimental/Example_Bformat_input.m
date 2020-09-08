%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Bformat input example
% ---------------------
% An example of how to render an input spherical harmonic (Ambisonic/
% B-Format) room impulse response (RIR) for a specific loudspeaker setup
% using HO-SIRR.
% 
% DEPENDENCES
%   Spherical-Harmonic-Transform Matlab library
%       https://github.com/polarch/Spherical-Harmonic-Transform
%   Higher-Order-Ambisonics Matlab library
%       https://github.com/polarch/Higher-Order-Ambisonics
%   Vector-Base-Amplitude-Panning
%       https://github.com/polarch/Vector-Base-Amplitude-Panning
%
% REFERENCES
%   [1] Favrot, S. and Buchholz, J.M., 2010. 
%       LoRA: A loudspeaker-based room auralization system. Acta Acustica 
%       united with Acustica, 96(2), pp.364-375.
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
%   Leo McCormack, 13/08/2019
%   leo.mccormack@aalto.fi 
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
clear all, close all, dbstop if error %#ok
addpath '..'
addpath '../_Simulated_Rooms_' '../_Stimuli_'

demo_order = 3;  
[input_stimulus, fs_in] = audioread('music__KickDrumClicky.wav');


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
% load('DTU_ls_dirs_deg.mat') 
% pars.ls_dirs_deg = ls_dirs_deg;
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
% tic, sirr_ls_rir = HOSIRR(sh_rir, pars); toc
% %audiowrite(['HOSIRR_o' num2str(demo_order) '.wav'], 0.9.*sirr_ls_rir, fs);
% 
% % Convolve the input stimulus with the rendered loudspeaker array RIR
% ir_length = size(sirr_ls_rir, 1);
% nLS = size(sirr_ls_rir, 2);  
% output = matrixConvolver(input_stimulus, reshape(sirr_ls_rir, [ir_length, 1, nLS]), size(input_stimulus,1));
% output = 0.99.*output./max(abs(output(:))); % normalise
% audiowrite(['BFormatDemo_HOSIRR_o' num2str(demo_order) '.wav'], output, fs, 'BitsPerSample', 24);

%% WIP: BINAURAL
% load HRIRs
pars.hrtf_sofa_path = '/Users/holdc1/Documents/data/HRTFs/Kemar_Aalto_2016/kemarhead_aalto2016.sofa';
[~, ~] = HOSIRR_bin(sh_rir, pars);

