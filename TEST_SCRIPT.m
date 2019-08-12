clear all, close all, dbstop if error %#ok
% this code is reliant on the following Matlab libraries:
%     https://github.com/polarch/Spherical-Harmonic-Transform
%     https://github.com/polarch/Higher-Order-Ambisonics
%
% Leo McCormack, 22/09/2018, leo.mccormack@aalto.fi

%% FIRST-ORDER SIRR EXAMPLE USING A-FORMAT MIC



%% HO-SIRR EXAMPLE USING EIGENMIKE (FREQUENCY-DEPENDENT ANALYSIS ORDER)
%(not implemented yet)




%% CHECKS
return;
%% HO-SIRR COMPARED WITH REFERENCE AND AMBISONICS
% A simulated Auditorium using LoRA, configured for a 36point t-design
% loudspeaker array [1].
LS_REFERENCE_FILE_NAME = 'Auditorium_36_Tdesign'; 

% USER PARAMETERS
pars.order = 3; 
pars.fs = 48e3; 
[~,dirs_rad] = getTdesign(8);
pars.ls_dirs_deg = dirs_rad*180/pi;    
pars.multires_winsize = [256];  
pars.multires_xovers = [];   
pars.RENDER_DIFFUSE = 1;
pars.BROADBAND_DIRECT = 1;   
pars.nBroadbandPeaks = 1;    
pars.decorrelationType = 'noise'; 
pars.maxDiffuseAnalysis_Hz = 6e3;  
pars.alpha_diff = 0.975;

% Encode reference audio into spherical harmonics
[refir, fs] = audioread([LS_REFERENCE_FILE_NAME '.wav']); % load input SHD IR, (ACN/N3D)
shir = refir * getRSH(pars.order, pars.ls_dirs_deg)';

% APPLY HO-SIRR
[sirr,sirr_ndiff,sirr_diff,pars,~] = HOSIRR(shir, pars, 0);
audiowrite([LS_REFERENCE_FILE_NAME '_HOSIRR_o' num2str(pars.order) '.wav'], 0.3.*sirr, fs);

% AMBISONICS
ambiIR = shir * getRSH(pars.order, pars.ls_dirs_deg)./size(pars.ls_dirs_deg, 1);

% Plots
subplot(1,3,1), s = surf(10*log10(abs(refir.')));
zlim([-50 0]), xlim([0 pars.fs/4]), xlabel('time'), ylabel('# channel'), title('ref')
s.EdgeColor = 'none';
colormap jet
subplot(1,3,2), s = surf(10*log10(abs(sirr.')));
zlim([-50 0]), xlim([0 pars.fs/4]), xlabel('time'), ylabel('# channel'), title('HO-SIRR')
s.EdgeColor = 'none';
colormap jet
subplot(1,3,3), s = surf(10*log10(abs(ambiIR.')));
zlim([-50 0]), xlim([0 pars.fs/4]), xlabel('time'), ylabel('# channel'), title('Ambisonics')
s.EdgeColor = 'none';
colormap jet
set(gcf,'Position',[40 40 1200 400])

clear pars;

% [1] Favrot, S. and Buchholz, J.M., 2010. LoRA: A loudspeaker-based room 
%     auralization system. Acta Acustica united with Acustica, 96(2), 
%     pp.364-375.


%% FIRST-ORDER SIRR CHECKS
pars.order = 1;  
pars.fs = 48e3; 
[~,dirs_rad] = getTdesign(2*(pars.order+1));
pars.ls_dirs_deg = dirs_rad*180/pi;
pars.multires_winsize = 256;  
pars.multires_xovers = [];   
pars.RENDER_DIFFUSE = 1;
pars.decorrelationType = 'noise'; 
pars.BROADBAND_DIRECT = 0;     
pars.maxDiffuseAnalysis_Hz = 6e3;  
pars.alpha_diff = 0.975;

% --- Single plane-wave input --- 
src_dir = [-45 -45];  % try adding more than 1 plane-wave, to see the first-order analysis break
shir = randn(pars.fs, size(src_dir,1)) * getRSH(pars.order, src_dir).';
[~,~,~,~,analysis] = HOSIRR(shir, pars, 0);

% In this case, we would expect:
% - [azimuth elevation] should correspond to 'src_dir'  
% - diffuseness should be 0
% - non-diffuse energy should be similar to input sound-field energy, and 
%   diffuse energy ~0
figure, subplot(4,1,1), imagesc(analysis.azim{1}*180/pi), colorbar, axis xy, caxis([-180 180]), title('azimuth (degrees)')
subplot(4,1,2), imagesc(analysis.elev{1}*180/pi), colorbar, axis xy, caxis([-90 90]), title('elevation (degrees)')
subplot(4,1,3), imagesc(10*log10(analysis.energy{1})), colorbar, axis xy, title('energy (dB)')
subplot(4,1,4), imagesc(analysis.diff{1}), colorbar, axis xy, caxis([0 1]), title('diffuseness')
figure, plot(10*log10(analysis.sf_energy{1})), hold on 
plot(10*log10(analysis.ndiff_energy{1})), hold on
plot(10*log10(analysis.diff_energy{1})) 
title('energy (dB)'), grid on, ylim([-40 20])
legend('sound-field', 'non-diffuse', 'diffuse')
figure, plot(analysis.diff{1}(1,:)), title('diffuseness'), grid on, ylim([0 1])

% --- Diffuse input --- 
[~, diff_dirs] = getTdesign(21); % approximate diffuse-field with 240 incoherent noise sources
shir = randn(pars.fs, size(diff_dirs,1)) * getRSH(pars.order, diff_dirs*180/pi).';
[~,~,~,~,analysis] = HOSIRR(shir, pars, 0);

% In this case, we would expect:
% - [azimuth elevation] should be random 
% - diffuseness should be close to 1
% - diffuse energy should be similar to the input sound-field energy, and
%   non-diffuse energy much lower than diffuse energy
figure, subplot(4,1,1), imagesc(analysis.azim{1}*180/pi), colorbar, axis xy, caxis([-180 180]), title('azimuth (degrees)')
subplot(4,1,2), imagesc(analysis.elev{1}*180/pi), colorbar, axis xy, caxis([-90 90]), title('elevation (degrees)')
subplot(4,1,3), imagesc(10*log10(analysis.energy{1})), colorbar, axis xy, title('energy (dB)')
subplot(4,1,4), imagesc(analysis.diff{1}), colorbar, axis xy, caxis([0 1]), title('diffuseness')
figure, plot(10*log10(analysis.sf_energy{1})), hold on 
plot(10*log10(analysis.ndiff_energy{1})), hold on
plot(10*log10(analysis.diff_energy{1})) 
title('energy (dB)'), grid on, ylim([-40 20])
legend('sound-field', 'non-diffuse', 'diffuse')
figure, plot(analysis.diff{1}(1,:)), title('diffuseness'), grid on, ylim([0 1])

clear pars


%% HO-SIRR CHECKS
pars.order = 3;  
pars.fs = 48e3; 
[~,dirs_rad] = getTdesign(2*(pars.order+1));
pars.ls_dirs_deg = dirs_rad*180/pi;
pars.multires_winsize = 256;  
pars.multires_xovers = [];   
pars.RENDER_DIFFUSE = 1;
pars.decorrelationType = 'noise'; 
pars.BROADBAND_DIRECT = 0;     
pars.maxDiffuseAnalysis_Hz = 6e3;  
pars.alpha_diff = 0.975;

% --- Three plane-wave input --- 
src_dir = [26 15; 153 -15; -116 -15]; % plane-waves landing in sectors: 1, 5, and 16
shir = randn(pars.fs, size(src_dir,1)) * getRSH(pars.order, src_dir).';
[~,~,~,~,analysis] = HOSIRR(shir, pars, 0);

% In this case, for sectors (1,5,16) we would expect:
% - [azimuth elevation] should correspond to 'src_dir' for their respective sectors 
% - These 3 sectors should have diffuseness values close to 0
for sectorIndex = [1,5,16]
    figure, subplot(4,1,1), imagesc(analysis.azim{1}(:,:,sectorIndex)*180/pi)
    colorbar, axis xy, caxis([-180 180]), title(['azimuth (degrees), sector: ' num2str(sectorIndex) ])
    subplot(4,1,2), imagesc(analysis.elev{1}(:,:,sectorIndex)*180/pi)
    colorbar, axis xy, caxis([-90 90]), title(['elevation (degrees), sector: ' num2str(sectorIndex) ])
    subplot(4,1,3), imagesc(10*log10(analysis.energy{1}(:,:,sectorIndex))), 
    colorbar, axis xy, title(['energy (dB), sector: ' num2str(sectorIndex) ])
    subplot(4,1,4), imagesc(analysis.diff{1}(:,:,sectorIndex))
    colorbar, axis xy, caxis([0 1]), title(['diffuseness, sector: ' num2str(sectorIndex) ]) 
end
% Whereas, sectors (2,20) we should expect:
% - [azimuth elevation] should be more random 
% - Higher diffuseness values
for sectorIndex = [2,20]
    figure, subplot(4,1,1), imagesc(analysis.azim{1}(:,:,sectorIndex)*180/pi)
    colorbar, axis xy, caxis([-180 180]), title(['azimuth (degrees), sector: ' num2str(sectorIndex) ])
    subplot(4,1,2), imagesc(analysis.elev{1}(:,:,sectorIndex)*180/pi)
    colorbar, axis xy, caxis([-90 90]), title(['elevation (degrees), sector: ' num2str(sectorIndex) ])
    subplot(4,1,3), imagesc(10*log10(analysis.energy{1}(:,:,sectorIndex))), 
    colorbar, axis xy, title(['energy (dB), sector: ' num2str(sectorIndex) ])
    subplot(4,1,4), imagesc(analysis.diff{1}(:,:,sectorIndex))
    colorbar, axis xy, caxis([0 1]), title(['diffuseness, sector: ' num2str(sectorIndex) ]) 
end 
% - Output non-diffuse energy should be similar to input sound-field 
%   energy, and also much higher than the diffuse energy
figure, plot(10*log10(analysis.sf_energy{1})), hold on 
plot(10*log10(analysis.ndiff_energy{1})), hold on
plot(10*log10(analysis.diff_energy{1})) 
title('energy (dB)'), grid on, ylim([-40 20])
legend('sound-field', 'non-diffuse', 'diffuse')

clear pars








