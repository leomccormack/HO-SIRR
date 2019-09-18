%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Energy Spectrogram Plotting Script 
% ----------------------------------
% A script for plotting energy spectrograms for different methods
% 
% DEPENDENCES
%     Spherical-Harmonic-Transform Matlab library
%         https://github.com/polarch/Spherical-Harmonic-Transform
%     Higher-Order-Ambisonics Matlab library
%         https://github.com/polarch/Higher-Order-Ambisonics
%     (Optional) SDM Toolbox
%         https://se.mathworks.com/matlabcentral/fileexchange/56663-sdm-toolbox
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
viewangs = [70 60];
display_freq_hz = 100;
enableSDMplots = exist('SDMbf', 'file');


%% Define and plot reference 
% this is a simulated Auditorium using LoRA [1] configured for DTU's AVIL 
% 64 channel anechoic listening room.
[ref, fs] = audioread('Auditorium_64ch_DTU_AVIL.wav');  
load('DTU_ls_dirs_deg.mat')
% [1] Favrot, S. and Buchholz, J.M., 2010. LoRA: A loudspeaker-based room 
%     auralization system. Acta Acustica united with Acustica, 96(2), 
%     pp.364-375.

% apply low-pass filter with specified cut-off, and compute the energy for
% each channel
[blp,alp]=butter(1, display_freq_hz/(fs/2)); % low-pass filter
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
pars.multires_winsize = 256;  
pars.multires_xovers = [];    
pars.RENDER_DIFFUSE = 1;
pars.BROADBAND_DIRECT = 1;   
pars.nBroadbandPeaks = 1;    
pars.decorrelationType = 'noise'; 
pars.maxDiffuseAnalysis_Hz = 6e3;  
pars.alpha_diff = 0.975;

% loop through the different input orders
for order = orders
    shir = ref * getRSH(order, ls_dirs_deg).'; 
    sirr = HOSIRR(shir, pars, 0); 

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


%% Mode-Matching Ambisonic decoder
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

%% SDM
if ~enableSDMplots; return; end 
shir = ref * getRSH(1, ls_dirs_deg).'; 

% convert shir to the ambisonic conventions used by SDM
shir = shir(:, [1 3 4 2]); % convert from ACN to WXYZ ordering
shir(:, 2:4) = shir(:, 2:4)./sqrt(3); % convert from N3D to SN3D normalisation

% configure and apply SDM, these values were taken from the demo scripts in
% the SDM toolbox. 
a = createSDMStruct('DefaultArray','Bformat','fs',fs,'winLen',15);
DOA = SDMbf(shir, a); 
P = shir(:,1);
ls_dirs_deg_w_r = [ls_dirs_deg ones(size(ls_dirs_deg,1),1)]; % add radius
s = createSynthesisStruct('lspLocs',ls_dirs_deg_w_r,'snfft',length(P) ,...
    'ShowArray',false,'fs',fs,'c',343,...
    'LFEchannel',[]); 
sdmir = synthesizeSDMCoeffs(P,DOA, s); 

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
title('SDM 1st order B-Format') 
set(gca,'fontsize', 14)
xlabel('Channel #')
ylabel('Time [ms]')
zlabel('Scaled energy [dB]')
set(gca,'Xtick',[size(ls_dirs_deg,1)/8:size(ls_dirs_deg,1)/8:size(ls_dirs_deg,1)]) %#ok
axis([1, size(ls_dirs_deg,1), 20, 400, -60, 0 ]) 
caxis([-60 -15])
colormap('jet')
%print -deps2c wfall_sdm.eps


%% SIRR 1st order no Diffuse Stream 
order = 1;

 % config
pars.chOrdering = 'ACN'; 
pars.normScheme = 'N3D';  
pars.ls_dirs_deg = ls_dirs_deg;
pars.fs = fs;  
pars.multires_winsize = 256;  
pars.multires_xovers = [];    
pars.RENDER_DIFFUSE = 0;  % disabled here
pars.BROADBAND_DIRECT = 1;   
pars.nBroadbandPeaks = 1;    
pars.decorrelationType = 'noise'; 
pars.maxDiffuseAnalysis_Hz = 6e3;  
pars.alpha_diff = 0.975;

shir = ref * getRSH(order, ls_dirs_deg).'; 
sirr = HOSIRR(shir, pars, 0); 

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
title(['SIRR no diffuse stream, order: ' num2str(order)]) 
set(gca,'fontsize', 14)
xlabel('Channel #')
ylabel('Time [ms]')
zlabel('Scaled energy [dB]')
set(gca,'Xtick',[size(ls_dirs_deg,1)/8:size(ls_dirs_deg,1)/8:size(ls_dirs_deg,1)]) %#ok
axis([1, size(ls_dirs_deg,1), 20, 400, -60, 0 ]) 
caxis([-60 -15])
colormap('jet')
%print -deps2c wfall_sirr.eps
