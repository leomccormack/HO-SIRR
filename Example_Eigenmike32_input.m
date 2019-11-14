%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Eigenmike32 input example
% -------------------------
% An example of how to render input Eigenmike32 room impulse responses for
% a specific loudspeaker arrangement using HO-SIRR.
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

demo_order = 4; 
nSH = int32((demo_order+1).^2);
[input_stimulus, fs_in] = audioread('music__KickDrumClicky.wav');

% Eigenmike32 angles, in [azimuth, elevation]:
mic_dirs_deg = ... 
    [0    21; 32   0; 0  -21; 328   0; 0    58; 45   35; 69    0; 45  -35; 0  -58; 
     315 -35; 291  0; 315 35; 91   69; 90   32; 90  -31; 89  -69; 180  21; 212  0; 
     180 -21; 148  0; 180 58; 225  35; 249   0; 225 -35; 180 -58; 135 -35; 111  0;
     135  35; 269 69; 270 32; 270 -32; 271 -69;];
mic_dirs_rad = mic_dirs_deg*pi/180;
R = 0.042; % the Eigenmike32 radius 
arrayType = 'rigid'; % the Eigenmike32 has a rigid baffle 


%% EITHER: CREATE DEMO EIGENMIKE ROOM IMPULSE RESPONSE (RIR)
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
audiowrite('Eigenmike32Demo_reference.wav', output, fs, 'BitsPerSample', 24);

% We then simulate an Eigenmike32, deriving the theoretical array transfer 
% functions for each loudspeaker direction in the DTU set-up
N_order = 30; % order of approximation 
Lfilt = 1024; % filter length
[h_mic, H_mic] = simulateSphArray(Lfilt, mic_dirs_rad, ls_dirs_rad, arrayType, R, N_order, fs);

% By convolving each loudspeaker RIR with the corresponding Eigenmike
% transfer functions, we essentially simulate an Eigenmike RIR 
% "measurement" of the reference room
mic_rir = matrixConvolver(ref_ls_rir, permute(h_mic, [1 3 2]), Lfilt/2);
mic_rir = 0.99.*mic_rir./max(abs(mic_rir(:)));

clear ls_dirs_deg ls_dirs_rad


%% OR: LOAD EIGENMIKE RIR FROM FILE 
%[mic_rir, fs] = audioread('PATH_TO_EIGENMIKE_RIR.wav');
disp(' * Replace "mic_rir" with your own Eigenmike room impulse response measurement')


%% EIGENMIKE RIR TO SPHERICAL HARMONIC RIR CONVERSION
Lfilt = 1024; % filter length
maxG_dB = 20; % maximum allowed gain amplification of the encoding filters
nMics = size(mic_dirs_rad,1);

% regularised inverse of the theoretical modal coefficients:
[H_filt, h_filt] = arraySHTfiltersTheory_radInverse(R, nMics, demo_order, Lfilt, fs, maxG_dB);
aziElev2aziPolar = @(dirs) [dirs(:,1) pi/2-dirs(:,2)]; % function to convert from azimuth-inclination to azimuth-elevation
Y_mics = sqrt(4*pi) * getSH(demo_order, aziElev2aziPolar(mic_dirs_rad), 'real'); % real SH matrix for microphone directions
h_filt_r = replicatePerOrder(h_filt/sqrt(4*pi), 2);
for kk=1:Lfilt/2+1
    M_mic2sh_radinv(:,:,kk) = diag(replicatePerOrder(H_filt(kk,:),2))*pinv(Y_mics); %#ok
end

% apply spherial harmonic transform
sh_rir_tmp = mic_rir * pinv(Y_mics).';
sh_rir = zeros(size(sh_rir_tmp,1)+Lfilt-1, nSH);
for i=1:nSH
    sh_rir(:,i) = fftconv(sh_rir_tmp(:,i), h_filt_r(:,i));  
end
sh_rir = 0.99.*sh_rir./max(abs(sh_rir(:)));
disp(' * Or, replace "sh_rir" with an encoding of your own Eigenmike room impulse response measurement')

 
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
audiowrite(['Eigenmike32Demo_HOSIRR_o' num2str(demo_order) '.wav'], output, fs, 'BitsPerSample', 24);

% The reference and HOSIRR renderings may then be compared over headphones using
% e.g. the free sparta_binauraliser VST plug-in [4] (which has a DTU setup preset).  
% Or, of course, you may also compare the two renders using the actual DTU 
% loudspeaker setup :-)


return; %Optional (and still experimental):
%% FIND USEABLE FREQUENCY RANGES FOR EACH SPHERICAL HARMONIC ORDER
% Obtain responses of the array for a dense grid of directions
[grid_azi, grid_elev] = meshgrid(-180:5:180, -85:5:85);
grid_dirs_deg = [grid_azi(:) grid_elev(:)];
grid_dirs_rad = grid_dirs_deg*pi/180;
Y_grid = sqrt(4*pi) * getSH(demo_order, aziElev2aziPolar(grid_dirs_rad), 'real'); % SH matrix for grid directions
N_order = 30; % order of approximation 
[~, H_array_sim] = simulateSphArray(Lfilt, mic_dirs_rad, grid_dirs_rad, arrayType, R, N_order, fs);

% Then compute the objective metrics, as described in [2,3]
evaluateSHTfilters(M_mic2sh_radinv, H_array_sim, fs, Y_grid);

% Note that the "Spatial correlation" is derived by comparing the patterns 
% of the array responses with the patterns of ideal spherical harmonics. 
% Here '1' means they are perfect, and '0' completely uncorrelated. 
% Therefore, the spatial aliasing frequency can be observed for each order,
% as the point where the spatial correlation tends towards 0. 
%
% The "Level difference" is then the  mean level difference over all 
% directions (diffuse level difference) between the ideal and simulated 
% components. One can observe that higher permitted amplification limits 
% [Max Gain (dB)] will result in noisier signals; however, this will also 
% result in a wider frequency range of useful spherical harmonic components
% at each order.

% After some soulsearching (and also a brief glance at these graphs), it is 
% perhaps reasonable to define the useable frequency ranges as
%    1st order up to 375 Hz, then:
%    2nd order   "   1e3 Hz,   "
%    3rd order   "   2e3 kHz,  "
%    4th order   "   6.5 kHz    "
%    3rd order   "   7.5 kHz   "
%    2nd order   "   8.5 kHz   "
%    1st order above 8.5 kHz   
pars.freqLimits = [375 1e3 2e3 6.5e3 7.5e3 8.5e3];
pars.procOrders = [1 2 3 4 3 2 1];

% Disclaimer: here we are using a theoretical description of the mic array. 
% Therefore, use these metrics just to get a rough idea of its performance. 
% Also note: strictly speaking, even the omni-directional (zeroth order) 
% component starts to alias above 9kHz. We only continue rendering past 
% this point for the sake of retaining the timbre. 


%% RENDER SH RIR TO TARGET LOUDSPEAKER LAYOUT
pars.alpha_diff = 0.975;

% same again, only... this time using frequency-dependent processing order:
tic, sirr_ls_rir = HOSIRR_orderPerBand(sh_rir, pars); toc 
%audiowrite(['HOSIRR_o' num2str(pars.maxOrder) '.wav'], 0.9.*sirr_ls_rir, fs);  

% Convolve the input stimulus with the rendered loudspeaker array RIR
ir_length = size(sirr_ls_rir, 1);
nLS = size(sirr_ls_rir, 2);  
output = matrixConvolver(input_stimulus, reshape(sirr_ls_rir, [ir_length, 1, nLS]), size(input_stimulus,1));
output = 0.99.*output./max(abs(output(:))); % normalise
audiowrite(['Eigenmike32Demo_HOSIRR_o' num2str(demo_order) 'orderPerBand.wav'], output, fs, 'BitsPerSample', 24);

 
