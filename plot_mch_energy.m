%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Energy Spectrogram Plotting Script 
% ----------------------------------
% A script for plotting energy spectrograms of loudspeaker array room 
% impulse responses (RIRs), using different rendering methods.
% 
% DEPENDENCES
%   Spherical-Harmonic-Transform Matlab library
%       https://github.com/polarch/Spherical-Harmonic-Transform
%   Higher-Order-Ambisonics Matlab library
%       https://github.com/polarch/Higher-Order-Ambisonics
%   Vector-Base-Amplitude-Panning
%       https://github.com/polarch/Vector-Base-Amplitude-Panning
%   Higher-Order-Ambisonics
%       https://github.com/polarch/Higher-Order-Ambisonics
%   (Optional) SDM Toolbox
%       https://se.mathworks.com/matlabcentral/fileexchange/56663-sdm-toolbox
%
% REFERENCES
%   [1] Favrot, S. and Buchholz, J.M., 2010. 
%       LoRA: A loudspeaker-based room auralization system. Acta Acustica 
%       united with Acustica, 96(2), pp.364-375.
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
%   Leo McCormack, 07/08/2019
%   leo.mccormack@aalto.fi 
%   Ville Pulkki, 23/05/2019
%   ville.pulkki@aalto.fi 
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
clear all, close all %#ok
addpath '_Simulated_Rooms_' '_Stimuli_'
viewangs = [70 60];
display_lpf_cutoff_freq_hz = 100;
enableSDMplots = exist('SDMbf', 'file');


%% PLOT REFERENCE
% This is a simulated Auditorium room using ODEON and LoRA [1], which is
% configured for DTU's AVIL anechoic listening room. 
% We will use this as a reference.
[ref, fs] = audioread('Auditorium_64ch_DTU_AVIL.wav');  
load('DTU_ls_dirs_deg.mat') 

% apply low-pass filter with the specified cut-off, and compute the energy for
% each channel
[blp,alp]=butter(1, display_lpf_cutoff_freq_hz/(fs/2)); % low-pass filter
e_ref=filter(blp,alp,ref.^2);
e_ref=10*log10(e_ref+eps);
maxx=max(max(e_ref));
e_ref=e_ref-maxx;
mean_e_ref = mean( mean(e_ref( (end-2048):end, :)) );
e_ref=max(-59.8,e_ref);

% plot resulting energy spectrogram
figure 
[m,n]=size(e_ref);
mesh([1:n]',[1:m]'/fs*1000,e_ref) %#ok
view(viewangs)
title('Reference') 
set(gca,'fontsize', 14)
xlabel('Channel #')
ylabel('Time [ms]')
zlabel('Scaled energy [dB]')
set(gca,'Xtick',[size(ls_dirs_deg,1)/8:size(ls_dirs_deg,1)/8:size(ls_dirs_deg,1)]) %#ok
axis([1, size(ls_dirs_deg,1), 20, 400, -60, 0 ])
caxis([-60 -15])
colormap('jet')
%print -deps2c wfall_ref.eps


%% HO-SIRR
orders = [1 3 5];

% config
pars.chOrdering = 'ACN'; 
pars.normScheme = 'N3D';  
pars.ls_dirs_deg = ls_dirs_deg;
pars.fs = fs;  
pars.multires_winsize = 128;  
pars.multires_xovers = [];  
%pars.multires_winsize = [512, 128, 64]; 
%pars.multires_xovers = [500, 2e3];  
pars.RENDER_DIFFUSE = 1;
pars.BROADBAND_FIRST_PEAK = 1;    
pars.BROADBAND_DIFFUSENESS = 1;
pars.maxDiffFreq_Hz = 3000;
pars.decorrelationType = 'noise';  
pars.alpha_diff = 0.5; 

% loop through the different input orders
for order = orders
    shir = ref * getRSH(order, ls_dirs_deg).'; 
    sirr = HOSIRR(shir, pars); %orderPerBand

    % apply same post-processing as with the reference
    e_sirr=filter(blp,alp,sirr.^2);
    e_sirr=10*log10(e_sirr+eps);
    maxx=max(max(e_sirr)); 
    e_sirr=e_sirr-maxx; 
    mean_e_sirr = mean( mean(e_sirr( (end-2048):end, :)) );

    % compensate for the change in maximum level, by normalising it to the
    % reference, based on the total energy in the last 2048 samples of the
    % response
    e_sirr = e_sirr-(mean_e_sirr-mean_e_ref);
    e_sirr=max(-59.8,e_sirr);

    % plot resulting energy spectrogram
    figure
    [m,n]=size(e_sirr);
    mesh([1:n]',[1:m]'/fs*1000,e_sirr) %#ok 
    view(viewangs)
    title(['SIRR order: ' num2str(order)]) 
    set(gca,'fontsize', 14)
    xlabel('Channel #')
    ylabel('Time [ms]')
    zlabel('Scaled energy [dB]')
    set(gca,'Xtick',[size(ls_dirs_deg,1)/8:size(ls_dirs_deg,1)/8:size(ls_dirs_deg,1)]) %#ok
    axis([1, size(ls_dirs_deg,1), 20, 400, -60, 0 ]) 
    caxis([-60 -15])
    colormap('jet')
    %print -deps2c wfall_sirr.eps
end

 
%% MODE-MATCHING AMBISONIC DECODER (MMD)
orders = [1 3 5];

% loop through the different input orders
for order = orders
    shir = ref * getRSH(order, ls_dirs_deg).'; 
    ambi = shir * pinv(getRSH(order, ls_dirs_deg).'); 

    % apply same post-processing as with the reference
    e_ambi=filter(blp,alp,ambi.^2);
    e_ambi=10*log10(e_ambi+eps);
    maxx=max(max(e_ambi)); 
    e_ambi=e_ambi-maxx; 
    mean_e_ambi = mean( mean(e_ambi( (end-2048):end, :)) );

    % compensate for the change in maximum level, by normalising it to the
    % reference, based on the total energy in the last 2048 samples of the
    % response
    e_ambi = e_ambi-(mean_e_ambi-mean_e_ref);
    e_ambi=max(-59.8,e_ambi);

    % plot resulting energy spectrogram
    figure
    [m,n]=size(e_ambi);
    mesh([1:n]',[1:m]'/fs*1000,e_ambi) %#ok 
    view(viewangs)
    title(['Ambisonics (MMD) order: ' num2str(order)]) 
    set(gca,'fontsize', 14)
    xlabel('Channel #')
    ylabel('Time [ms]')
    zlabel('Scaled energy [dB]')
    set(gca,'Xtick',[size(ls_dirs_deg,1)/8:size(ls_dirs_deg,1)/8:size(ls_dirs_deg,1)]) %#ok
    axis([1, size(ls_dirs_deg,1), 20, 400, -60, 0 ]) 
    caxis([-60 -15])
    colormap('jet')
    %print -deps2c wfall_ambi.eps
end


%% SPATIAL DECOMPOSITION METHOD (B-FORMAT VARIANT)
if ~enableSDMplots; return; end 
shir = ref * getRSH(1, ls_dirs_deg).'; 

% convert shir to the ambisonic conventions used by SDM
shir = shir(:, [1 4 2 3]); % convert from ACN to WXYZ ordering
shir(:, 2:4) = shir(:, 2:4)./sqrt(3); % convert from N3D to SN3D normalisation

% configure and apply SDM 
% (note that these values were taken from the b-format demo script in the 
% SDM toolbox) - Except a window length of 1, instead of 15
a = createSDMStruct('DefaultArray','Bformat','fs',fs,'winLen',1);
DOA = SDMbf(shir, a); 
P = shir(:,1);
ls_dirs_deg_w_r = [ls_dirs_deg ones(size(ls_dirs_deg,1),1)]; % add radius
s = createSynthesisStruct('lspLocs',ls_dirs_deg_w_r,'snfft',length(P) ,...
    'ShowArray',false,'fs',fs,'c',343,...
    'LFEchannel',[]); 
sdmir = synthesizeSDMCoeffs(P,DOA, s); 

% figure, subplot(2,1,1)
% plot(shir(1:end/8,1)), title('input pressure')
% subplot(2,1,2)
% plot(sum(sdmir(1:end/8,:),2)), title('combined response')

% apply same post-processing as with the reference
e_sdmir=filter(blp,alp,sdmir.^2);
e_sdmir=10*log10(e_sdmir+eps);
maxx=max(max(e_sdmir)); 
e_sdmir=e_sdmir-maxx; 
mean_e_sirr = mean( mean(e_sdmir( (end-2048):end, :)) );

% compensate for the change in maximum level, by normalising it to the
% reference, based on the total energy in the last 2048 samples of the
% response
e_sdmir = e_sdmir-(mean_e_sirr-mean_e_ref);
e_sdmir=max(-59.8,e_sdmir);

% plot resulting energy spectrogram
figure
[m,n]=size(e_sdmir);
mesh([1:n]',[1:m]'/fs*1000,e_sdmir) %#ok 
view(viewangs)
title('SDM (B-Format) order: 1') 
set(gca,'fontsize', 14)
xlabel('Channel #')
ylabel('Time [ms]')
zlabel('Scaled energy [dB]')
set(gca,'Xtick',[size(ls_dirs_deg,1)/8:size(ls_dirs_deg,1)/8:size(ls_dirs_deg,1)]) %#ok
axis([1, size(ls_dirs_deg,1), 20, 400, -60, 0 ]) 
caxis([-60 -15])
colormap('jet')
%print -deps2c wfall_sdm.eps


%% SIRR 1ST ORDER NO DIFFUSE STREAM 
order = 1;

 % config
pars.chOrdering = 'ACN'; 
pars.normScheme = 'N3D';  
pars.ls_dirs_deg = ls_dirs_deg;
pars.fs = fs;  
pars.multires_winsize = 16; % more in line with SDM
pars.multires_xovers = [];    
pars.RENDER_DIFFUSE = 0;  % disabled here
pars.BROADBAND_FIRST_PEAK = 1;    
pars.BROADBAND_DIFFUSENESS = 1;
pars.maxDiffFreq_Hz = 3000; 
pars.decorrelationType = 'noise';  
pars.alpha_diff = 0.5; 

shir = ref * getRSH(order, ls_dirs_deg).'; 
sirr = HOSIRR(shir, pars); 

% apply same post-processing as with the reference
e_sirr=filter(blp,alp,sirr.^2);
e_sirr=10*log10(e_sirr+eps);
maxx=max(max(e_sirr)); 
e_sirr=e_sirr-maxx; 
mean_e_sirr = mean( mean(e_sirr( (end-2048):end, :)) );

% compensate for the change in maximum level, by normalising it to the
% reference, based on the total energy in the last 2048 samples of the
% response
e_sirr = e_sirr-(mean_e_sirr-mean_e_ref);
e_sirr=max(-59.8,e_sirr);

% plot resulting energy spectrogram
figure 
[m,n]=size(e_sirr);
mesh([1:n]',[1:m]'/fs*1000,e_sirr) %#ok 
view(viewangs)
title(['SIRR without diffuse stream order: ' num2str(order)]) 
set(gca,'fontsize', 14)
xlabel('Channel #')
ylabel('Time [ms]')
zlabel('Scaled energy [dB]')
set(gca,'Xtick',[size(ls_dirs_deg,1)/8:size(ls_dirs_deg,1)/8:size(ls_dirs_deg,1)]) %#ok
axis([1, size(ls_dirs_deg,1), 20, 400, -60, 0 ]) 
caxis([-60 -15])
colormap('jet')
%print -deps2c wfall_sirr.eps


%% HO-SIRR (ORDER PER BAND VERSION)
order = 4;
pars.freqLimits = [375 1e3 2e3 6.5e3 7.5e3 8.5e3];
pars.procOrders = [1 2 3 4 3 2 1];

% config
pars.chOrdering = 'ACN'; 
pars.normScheme = 'N3D';  
pars.ls_dirs_deg = ls_dirs_deg;
pars.fs = fs;  
pars.multires_winsize = 128;  
pars.multires_xovers = [];  
%pars.multires_winsize = [512, 128, 64]; 
%pars.multires_xovers = [500, 2e3];  
pars.RENDER_DIFFUSE = 1;
pars.BROADBAND_FIRST_PEAK = 1;     
pars.decorrelationType = 'noise';  
pars.alpha_diff = 0.975; 

shir = ref * getRSH(order, ls_dirs_deg).'; 
sirr = HOSIRR_orderPerBand(shir, pars); 

% apply same post-processing as with the reference
e_sirr=filter(blp,alp,sirr.^2);
e_sirr=10*log10(e_sirr+eps);
maxx=max(max(e_sirr)); 
e_sirr=e_sirr-maxx; 
mean_e_sirr = mean( mean(e_sirr( (end-2048):end, :)) );

% compensate for the change in maximum level, by normalising it to the
% reference, based on the total energy in the last 2048 samples of the
% response
e_sirr = e_sirr-(mean_e_sirr-mean_e_ref);
e_sirr=max(-59.8,e_sirr);

% plot resulting energy spectrogram
figure
[m,n]=size(e_sirr);
mesh([1:n]',[1:m]'/fs*1000,e_sirr) %#ok 
view(viewangs)
title('SIRR step-wise order (Eigenmike ranges)') 
set(gca,'fontsize', 14)
xlabel('Channel #')
ylabel('Time [ms]')
zlabel('Scaled energy [dB]')
set(gca,'Xtick',[size(ls_dirs_deg,1)/8:size(ls_dirs_deg,1)/8:size(ls_dirs_deg,1)]) %#ok
axis([1, size(ls_dirs_deg,1), 20, 400, -60, 0 ]) 
caxis([-60 -15])
colormap('jet')
%print -deps2c wfall_sirr.eps


 