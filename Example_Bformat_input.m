%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Bformat input example
% ---------------------
% 
% DEPENDENCES
%     Spherical-Harmonic-Transform Matlab library
%         https://github.com/polarch/Spherical-Harmonic-Transform
%     Higher-Order-Ambisonics Matlab library
%         https://github.com/polarch/Higher-Order-Ambisonics
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
%   Leo McCormack, 13/08/2019
%   leo.mccormack@aalto.fi 
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
clear all, close all, dbstop if error %#ok


%% Define Input IR
% For demonstration purposes, this is a simulated Auditorium using LoRA [1]
% configured for a 36point t-design array [1].
demo_order = 3; % order to encode Auditorium
[refir, fs] = audioread('Auditorium_36_Tdesign.wav');  
[~,t_dirs_rad] = getTdesign(8); % directions at which the Auditorium was simulated
t_dirs_deg = t_dirs_rad*180/pi; 
shir = refir * getRSH(demo_order, t_dirs_deg)'; % encode to spherical harmonics
% [1] Favrot, S. and Buchholz, J.M., 2010. LoRA: A loudspeaker-based room 
%     auralization system. Acta Acustica united with Acustica, 96(2), 
%     pp.364-375.

% !!!
% Replace 'shir' with your own B-format/Ambisonics/spherical harmonic signals
% !!!


%% User parameters
% Specify input signal channel ordering and normalisation conventions
% Note: 
%     Ambix -> ACN/SN3D
%     FuMa  -> WXYZ/SN3D (first-order only. Note this is also without the 
%     1/sqrt(2) scaling on the omni component)
pars.chOrdering = 'ACN'; % 'ACN', or 'WXYZ' (deprecated/first-order only)
pars.normScheme = 'N3D'; % 'N3D', or 'SN3D'
% Define your loudspeaker directions here, [azim elev] convention, in 
% DEGREES; nLoudspeakers x 2
pars.ls_dirs_deg = [ 0 0; 45 0; -45 0; 90 0; -90 0; 135 0; -135 0; 45 35; -45 35; 90 35; -90 35; 0 90;]; 


%% Render Input to Loudspeakers using HO-SIRR
% Optionally, the below parameters may also be changed.
pars.fs = fs;  
pars.multires_winsize = 256; % windowing size, in samples, note HOSIRR employs 50% overlap
pars.multires_xovers = [];   
% or e.g.:
% 1024 hop size < 500Hz < 128 hop size < 2kHz < 32 hop size
%pars.multires_winsize = [1024, 128, 32]; 
%pars.multires_xovers = [500, 2e3];   
pars.RENDER_DIFFUSE = 1;
pars.BROADBAND_DIRECT = 1;   
pars.nBroadbandPeaks = 1;    
pars.decorrelationType = 'noise'; 
pars.maxDiffuseAnalysis_Hz = 6e3;  
pars.alpha_diff = 0.975;
 
% Render
[sirr,~,~,pars,~] = HOSIRR(shir, pars, 0);
audiowrite(['demo_render_HOSIRR_o' num2str(pars.maxOrder) '.wav'], 0.3.*sirr, fs);



