%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Upmixing B-format example
% -------------------------
% An example of how to upmix a LOWER-order ambisonic room impulse response 
% (RIR) to a HIGHER-order ambisonic RIR using HO-SIRR.
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
addpath '_Simulated_Rooms_'  

input_order = 1;
output_order = 3;  
  

%% EITHER: CREATE DEMO INPUT RIR
% For demonstration purposes, this is a simulated Auditorium using LoRA [1]
% configured for DTU's AVIL anechoic listening room.
[ref_ls_rir, fs] = audioread('Auditorium_64ch_DTU_AVIL.wav');  
load('DTU_ls_dirs_deg.mat')

% encode into the spherical harmonic domain at "input_order"
nLS = size(ls_dirs_deg,1);
input_sh_rir = ref_ls_rir * sqrt(4*pi/nLS) * getRSH(input_order, ls_dirs_deg)'; 
input_sh_rir = 0.99.*input_sh_rir./max(abs(input_sh_rir(:)));
audiowrite(['UpmixingDemo_HOSIRR_input_o' num2str(input_order) '.wav'], input_sh_rir, fs, 'BitsPerSample', 24);
clear ls_dirs_deg;


%% OR: LOAD ROOM IMPULSE RESPONSE FROM FILE 
%[input_sh_rir, fs] = audioread('PATH_TO_SH_RIR.wav');
disp(' * Replace "input_sh_rir" with your own Ambisonic/B-Format room impulse response measurement')


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


%% RENDER SH RIR TO T-DESIGN AND ENCODE BACK INTO SPHERICAL HARMONIC DOMAIN
% suitable t-design for this target order
[~,ls_dirs_rad] = getTdesign(2*output_order);
pars.ls_dirs_deg = ls_dirs_rad * 180/pi;

% render
tic, sirr_ls_rir = HOSIRR(input_sh_rir, pars); toc

% encode
nLS = size(pars.ls_dirs_deg,1);
output_sh_rir = sirr_ls_rir * sqrt(4*pi/nLS) * getRSH(output_order, pars.ls_dirs_deg)'; 
output_sh_rir = 0.99.*output_sh_rir./max(abs(output_sh_rir(:)));
audiowrite(['UpmixingDemo_HOSIRR_input_o' num2str(input_order) '_output_o' num2str(output_order) '.wav'], output_sh_rir, fs, 'BitsPerSample', 24);
 
