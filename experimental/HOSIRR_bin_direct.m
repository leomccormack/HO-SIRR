function [lsir, lsir_ndiff, lsir_diff, pars, analysis] = HOSIRR_bin(shir, pars)
% Higher-Order Spatial Impulse Response Rendering (HOSIRR)
% --------------------------------------------------------
% Multi-resolution Higher-order Spatial Impulse Response Rendering
% 
% DEPENDENCES
%   Spherical-Harmonic-Transform Matlab library
%       https://github.com/polarch/Spherical-Harmonic-Transform
%   Higher-Order-Ambisonics Matlab library
%       https://github.com/polarch/Higher-Order-Ambisonics
%
% INPUT ARGUMENTS
%   shir                       : spherical harmonic domain impulse response
%                                [signalLength  x (order+1)^2]
%   pars.chOrdering            : {'ACN','WXYZ'} channel ordering convention
%                                Note: 'WXYZ' refers to "FuMa", which is 
%                                first-order only.
%   pars.normScheme            : {'N3D','SN3D'}, normalisation convention. 
%                                fully normalised (N3D), or semi-normalised
%                                (SN3D) conventions are supported.
%   pars.fs                    : sample rate
%   pars.ls_dirs_deg           : loudspeaker array directions in DEGREES 
%                                [azi elev] convention
%   pars.multires_winsize      : [winsize_low winsize_2 ... winsize_high]
%   pars.multires_xovers       : [xover_low ... xover_high] 
%   pars.RENDER_DIFFUSE        : {0,1} 0 off, 1 new diffuse stream via
%                                ambisonic decoding of diffuseness scaled
%                                sector signals
%   pars.BROADBAND_FIRST_PEAK  : {0,1} 0 off, 1 broadband analysis for 
%                                direct
%   pars.BROADBAND_DIFFUSENESS : {0,1} 0 diffuseness computed per
%                                frequency bin, 1 diffuseness computed up
%                                to "maxDiffFreq_Hz", and replicated for
%                                all bins
%   pars.NORMALIZE_INPUT       : Normalize input to max(|insig|) = 1
%   pars.maxDiffFreq_Hz        : frequency up to which to compute the
%                                diffuseness parameter for
%   pars.alpha_diff            : one-pole alpha value for smoothing diff 
%                                y(n) = alpha*y(n-1) + (1-alpha)*x(n)
%   pars.decorrelationType     : {'phase','noise'}  
%
% OUTPUT ARGUMENTS
%   lsir                       : loudspeaker impulse responses    
%   lsir_ndiff                 : non-diffuse stream only
%   lsir_diff                  : diffuse stream only
%   pars                       : parameter struct used for the processing
%   analysis                   : {nRes}[nBins x nFrames x nSectors],
%                                historic DoA estimates (.azim, .elev), 
%                                .energy and diffuseness (.diff) estimates
%
% REFERENCES
%   [1] McCormack, L., Pulkki, V., Politis, A., Scheuregger, O. and Marschall,
%       M., (2020). "Higher-Order Spatial Impulse Response Rendering:
%       Investigating the Perceived Effects of Spherical Order, Dedicated
%       Diffuse Rendering, and Frequency Resolution". Journal of the Audio
%       Engineering Society, 68(5), pp.338-354.
%   [2] McCormack, L., Politis, A., Scheuregger, O., and Pulkki, V. 2019.
%       "Higher-order processing of spatial impulse responses". In 
%       Proceedings of the 23rd International Congress on Acoustics, 
%       9--13 September 2019 in Aachen, Germany.
%   [3] Politis, A. and Pulkki, V., 2016. "Acoustic intensity, energy-
%       density and diffuseness estimation in a directionally-constrained 
%       region". arXiv preprint arXiv:1609.03409.
%   [4] Merimaa, J. and Pulkki, V., 2005. "Spatial impulse response
%       rendering I: Analysis and synthesis". Journal of the Audio 
%       Engineering Society, 53(12), pp.1115-1127.
%   [5] Pulkki, V. and Merimaa, J., 2006. "Spatial impulse response
%       rendering II: Reproduction of diffuse sound and listening tests". 
%       Journal of the Audio Engineering Society, 54(1/2), pp.3-20. 
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
%   Christoph Hold, 2022
%   Christoph.Hold@aalto.fi
%   Leo McCormack, 22/09/2018
%   leo.mccormack@aalto.fi
%   Archontis Politis, 12/06/2018
%   archontis.politis@aalto.fi
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

nSH = size(shir,2);
pars.order = sqrt(nSH)-1;

%%% Defaults/Warnings/Errors: 
if ( (pars.order>1) && strcmp(pars.chOrdering, 'WXYZ') )
    disp('Warning: WXYZ/FuMa is first-order only. Input signals have been truncated to first-order')
    pars.maxOrder = 1;
    nSH = int32((pars.maxOrder+1).^2);
    shir = shir(:,1:nSH);
end
if (~isfield(pars, 'fs')), error('Please specify "fs"'); end
%if (~isfield(pars, 'ls_dirs_deg')) 
%    error('Please specify "ls_dirs_deg", in degrees')
%end 
if (~isfield(pars, 'multires_winsize')), pars.multires_winsize = 512; end
if (~isfield(pars,'multires_xovers') || isempty(pars.multires_xovers))
    nRes = 1;
elseif (length(pars.multires_winsize)~=length(pars.multires_xovers)+1)
    error('The number of specified window sizes does not match the number of resolutions')
else
    nRes = length(pars.multires_winsize);
end 


%%% HO-SIRR
disp('HOSIRR Configuration:'), pars %#ok

% convert to 'ACN/N3D' Ambisonic convention
if strcmp(pars.chOrdering, 'WXYZ')
    shir = convert_N3D_Bformat(shir, 'b2n'); 
elseif (strcmp(pars.chOrdering, 'ACN') && strcmp(pars.normScheme, 'SN3D'))
    shir = convert_N3D_SN3D(shir, 'sn2n');
end

lSig = size(shir,1);

if isfield(pars,'NORMALIZE_INPUT') && pars.NORMALIZE_INPUT
    % normalise input to max(|insig|) = 1
    shir = shir./(max(abs(shir(:,1))));
end

% extract the first highest peak
if pars.BROADBAND_FIRST_PEAK
    shir_tmp = shir; 
    
    % find the index of highest peak in the omni
    [~, peak_ind] = max(abs(shir_tmp(:,1)).^2);

    % calculate window 
    dirwinsize = 64;
    direct_win = zeros(lSig,1);
    direct_win(peak_ind-dirwinsize/2:peak_ind+dirwinsize/2,1) = hanning(dirwinsize+1);

    % extract peak from shir
    shir_direct = repmat(direct_win, [1 nSH]).*shir;  
    shir_tmp = repmat(1-direct_win, [1 nSH]).*shir_tmp;  

    shir = shir_tmp; % remove peak before running the main loop
    clear shir_tmp
end

%%% Intialisations
% zero pad the signal's start and end for STFT
maxWinsize = max(pars.multires_winsize);
if maxWinsize == 1
    shir_pad = shir;
else 
    shir_pad = [zeros(maxWinsize/2, nSH); shir; zeros(maxWinsize*2, nSH)];
end

%% INTIT HRTFS
if (~isfield(pars, 'hrirs'))
pars = prepareHRIRs(pars, pars.hrtf_sofa_path);
end

hrirs = pars.hrirs;
lenHrirs = pars.lenHrirs;


pars.hrtf_vbapTableRes = [2 5]; % resolution azim/elev in degs
vbapTable = getGainTable(pars.hrir_dirs_deg, pars.hrtf_vbapTableRes);
numTable = size(vbapTable,1);

fprintf('\nComputing interpolation gain table ...... ')
pars.hrtf_vbapTableComp = zeros(numTable, 3); 
pars.hrtf_vbapTableIdx = ones(numTable, 3); % ones! not zeros!
% compress table by keeping only the non-zero gains and their indices
for nt=1:numTable
    temp_nt = vbapTable(nt,:);
    gains_nt = temp_nt(temp_nt>0);
    idx_nt = find(temp_nt>0); 
    pars.hrtf_vbapTableComp(nt,1:length(gains_nt)) = gains_nt;
    pars.hrtf_vbapTableIdx(nt,1:length(idx_nt)) = idx_nt;
end
clear vbapTable
pars.hrtf_vbapTableComp = pars.hrtf_vbapTableComp ./ (sum(pars.hrtf_vbapTableComp,2)*ones(1,3));


% Loudspeaker locations and VBAP gain table
%ls_dirs_deg = pars.ls_dirs_deg;
%nLS = size(ls_dirs_deg,1);
%vbap_gtable_res = [2 2]; % azi, elev, step sizes in degrees
%gtable = getGainTable(ls_dirs_deg, vbap_gtable_res); 


% Sector design
if (~isfield(pars, 'pattern')), pars.pattern = 'pwd'; end
if (~isfield(pars, 'ENABLE_DIFF_PR')), pars.ENABLE_DIFF_PR = true; end

[~,sec_dirs_rad] = getTdesign(2*(max([1, pars.order-1])));  % not consistent with HOSIRR

% Rotation for alignment
secPos = unitSph2cart(sec_dirs_rad);
R = calculateRotationMatrix(secPos(1, :), [1, 0, 0]);
secPos = secPos * R;  % rotate to sec0 = [0,0]
[secAzi, secEle, secR] = cart2sph(secPos(:,1), secPos(:,2), secPos(:,3));
sec_dirs_rad = [secAzi, secEle];

A_xyz = computeVelCoeffsMtx(pars.order-1);
[pars.sectorCoeffs, pars.secNorms, sec_dirs_rad] = ...
    computeSectorCoeffs(pars.order-1, A_xyz, pars.pattern, [], sec_dirs_rad);
% amplitude normalisation term
beta_A = pars.secNorms(1);
% energy normalisation term
beta_E = pars.secNorms(2);
% Perfect Reconstruction of diffuse component
if pars.ENABLE_DIFF_PR
    switch pars.pattern
        case 'cardioid'
            b_n_an = beamWeightsCardioid2Spherical(pars.order-1);
        case 'maxRE'
            b_n_an = beamWeightsMaxEV(pars.order-1);
        case 'pwd'
            b_n_an = beamWeightsHypercardioid2Spherical(pars.order-1);
    end
    c_n = b_n_an ./ sqrt((2*(0:pars.order-1).'+1) / (4*pi));
    c_n_res = 1./(c_n./c_n(1));
    % TODO
    % Bandlimit
else
    c_n_res = ones(pars.order, 1);
end


if pars.order~=1
    pars.sectorDirs = sec_dirs_rad;
end 
numSec = size(pars.sectorCoeffs,2)/4;

% divide signal to frequency regions
if (nRes>1)
    xover_order = 500;
    shir_res = divide2regions(shir_pad, pars.fs, pars.multires_xovers, xover_order);
else
    xover_order = 0;
    shir_res = shir_pad;
end
clear insig_pad
  
% time-frequency processing for each frequency region
assert(maxWinsize <= lenHrirs)  % for now
maxfftsize = 2*lenHrirs;
lsir_res_ndiff = zeros(lSig + 2*maxfftsize + xover_order, 2, nRes);  % 2 ears
lsir_res_diff = zeros(lSig + 2*maxfftsize + xover_order, 2, nRes); 

%%% Multi-resolution processing
for nr = 1:nRes 
    disp(['Processing frequency region no. ' num2str(nr)]);
    winsize = pars.multires_winsize(nr);
    hopsize = winsize/2; % half the window size time-resolution
    %nBins_anl = winsize/2+1; % nBins used for analysis
    nBins_anl = 1*winsize/2+1;  % Zero pad / interpolate input with this
    %disp(nBins_anl)
    % try to prevent something stupid
    assert(nBins_anl >= winsize/2+1)
    assert(nBins_anl <= lenHrirs/2+1)
    
    % Assumes win <= hrir
    fftsize_syn = 2*lenHrirs; % double the window size for FD convolution, TODO: more than necessary.
    nBins_syn = fftsize_syn/2 + 1; % nBins used for synthesis 
    pars.centerfreqs_anl = (0: nBins_anl -1)'*pars.fs/((nBins_anl-1)*2);
    pars.hrtf_centerfreqs = (0: lenHrirs/2+1 -1)'*pars.fs/lenHrirs;
    if pars.BROADBAND_DIFFUSENESS
        if nr==nRes
            [~,maxDiffFreq_Ind] = min(abs(pars.centerfreqs_anl-min(pars.maxDiffFreq_Hz)));
        else
            [~,maxDiffFreq_Ind] = min(abs(pars.centerfreqs_anl-min(pars.maxDiffFreq_Hz,pars.multires_xovers(nr))));  
        end 
    end   
    
    % Prepare HRTFs
    hrtfs = fft(hrirs, lenHrirs, 1);
    %hrtfs_syn = fft(hrir, fftsize, 1);
    if mod(size(hrirs, 1), 2)
        error('not implemented')
    else  % even
        hrtfs = hrtfs(1:lenHrirs/2+1, :, :);  % pos half
    end
    %pars.hrtfs = hrtfs;
    pars.hrtf_mag = abs(hrtfs);

    [x_ls, y_ls, z_ls] = sph2cart( pars.sectorDirs(:,1)*pi/180,...
                               pars.sectorDirs(:,2)*pi/180, 1);
    [x_hrfts, y_hrtfs, z_hrtfs] = sph2cart(pars.hrir_dirs_deg(:, 1)*pi/180,...
                                           pars.hrir_dirs_deg(:, 2)*pi/180, 1);
    ls_proj = [x_hrfts, y_hrtfs, z_hrtfs] * [x_ls, y_ls, z_ls].';
    [~, d_min_k] = max(ls_proj);
    hrtfs_sec = hrtfs(:,:,d_min_k);
     
    % storage for estimated parameters
    analysis.azim{nr} = nan(nBins_anl, ceil(lSig/hopsize), numSec);
    analysis.elev{nr} = nan(nBins_anl, ceil(lSig/hopsize), numSec);
    analysis.energy{nr} = nan(nBins_anl, ceil(lSig/hopsize), numSec);
    analysis.diff{nr} = nan(nBins_anl, ceil(lSig/hopsize), numSec);
     
    % extended energy analysis
    analysis.sf_energy{nr} = nan(ceil(lSig/hopsize),1);
    analysis.ndiff_energy{nr} = nan(ceil(lSig/hopsize),1);
    analysis.diff_energy{nr} = nan(ceil(lSig/hopsize),1);
    
    % transform window (hanning)
    x = 0:(winsize-1);
    win = sin(x.*(pi/winsize))'.^2;
      
    % diffuse stream rendering intialisations
    switch pars.RENDER_DIFFUSE
        case 0
            % No diffuse stream
        case 1
            % New SIRR diffuse stream, based on scaling the sector signals 
            % with diffuseness estimates and then re-encoding them into 
            % SHs, and then decoding them binaurally.
            if pars.order==1
                %D_ls = sqrt(4*pi/nLS).*getRSH(pars.order, ls_dirs_deg).';
                D_bin = getAmbisonic2BinauralFilters_magls_hfcont(...
                    permute(hrtfs, [2 3 1]), pars.hrir_dirs_deg, pars.order, pars.fs, 1500, pars.hrirs_weights);
            else 
                Y_enc = diag(replicatePerOrder(c_n_res)) * ...
                    getRSH(pars.order-1, pars.sectorDirs*180/pi); % encoder
                %D_ls = sqrt(4*pi/nLS).*getRSH(pars.order-1, ls_dirs_deg).';   
                D_bin = getAmbisonic2BinauralFilters_magls_hfcont(...
                    permute(hrtfs, [2 3 1]), pars.hrir_dirs_deg, pars.order-1, pars.fs, 1500, pars.hrirs_weights);
            end
            D_bin_interp = interpolateFilters(D_bin, fftsize_syn);
    end 
     
    % diffuseness averaging buffers 
    if ~pars.BROADBAND_DIFFUSENESS
        prev_intensity = zeros(nBins_anl,3,numSec);
        prev_energy = zeros(nBins_anl,numSec);
    else 
        prev_intensity = zeros(3,numSec);
        prev_energy = zeros(numSec,1);
    end
      
    % analysed parameters for each sector
    azim = zeros(nBins_anl,numSec);
    elev = zeros(nBins_anl,numSec);
    diffs = zeros(nBins_anl,numSec); 
    
    %%% Main processing loop
    idx = 1;
    framecount = 1;
    progress = 1;
    nFrames = ceil((lSig + maxWinsize)/hopsize)+1;
    
    while idx + maxWinsize <= lSig + 2*maxWinsize  
        % Window input and transform to frequency domain
        insig_win = win*ones(1,nSH) .* shir_res(idx+(0:winsize-1),:,nr);
        %inspec = fft(insig_win);
        %inspec = inspec(1:nBins_anl,:); % keep up to nyquist
        inspec = fft(insig_win, (nBins_anl-1)*2);  % interpolates if needed
        inspec = inspec(1:nBins_anl,:); % keep up to nyquist
         
        %%% SIRR ANALYSIS %%%
        outspec_ndiff = 0;
        outspec_diff = 0;
        W_S = pars.sectorCoeffs;   
        s_ana = inspec*W_S; 
        for n=1:numSec 
            % weighted pressure-velocity signals for this sector
            WXYZ_sec = s_ana(:,4*(n-1) + (1:4));   
            
            % Compute Intensity vector for each frequency bin
            I = real(conj(WXYZ_sec(:,1)*ones(1,3)) .* WXYZ_sec(:,2:4));  
            [azim(:,n), elev(:,n)] = cart2sph(I(:,1), I(:,2), I(:,3));  
                 
            if pars.BROADBAND_DIFFUSENESS
                % Compute broad-band active-intensity vector
                pvCOV = (WXYZ_sec(1:maxDiffFreq_Ind,:)'*WXYZ_sec(1:maxDiffFreq_Ind,:)); 
                I_diff = real(pvCOV(2:4,1));
                energy = 0.5.*real(trace(pvCOV));  % cast back to real

                % Estimating and time averaging of boadband diffuseness
                diff_intensity = (1-pars.alpha_diff).*I_diff + pars.alpha_diff.*prev_intensity(:,n);
                diff_energy = (1-pars.alpha_diff).*energy + pars.alpha_diff.*prev_energy(n); 
                prev_intensity(:,n) = diff_intensity;
                prev_energy(n) = diff_energy; 
                diffs(:,n) = 1 - sqrt(sum(diff_intensity.^2)) ./ (diff_energy + eps); 
            else  
                energy = 0.5.*sum(abs(WXYZ_sec).^2,2); 
            
                % Time averaging of intensity-vector for the diffuseness
                % estimate per bin
                diff_intensity = (1-pars.alpha_diff).*I + pars.alpha_diff.*prev_intensity(:,:,n);
                diff_energy = (1-pars.alpha_diff).*energy + pars.alpha_diff.*prev_energy(:,n); 
                diffs(:,n) = 1 - sqrt(sum(diff_intensity.^2,2)) ./ (diff_energy + eps); 
                prev_intensity(:,:,n) = diff_intensity;
                prev_energy(:,n) = diff_energy;
                % freq smoothing  (avoid sharp modulation)
                if true
                    diffs(:,n) = medfilt1(diffs(:,n), 4, [], 1, 'truncate');
                    %diffs(:,n) = sgolayfilt(diffs(:,n), 3, 5, [], 1);
                end
                %assert(all(diffs(:,n)<=1.001))
                %assert(all(diffs(:,n)>=0))
            end


            
            % storage for estimated parameters over time
            analysis.azim{nr}(:,framecount,n) = azim(:,n);
            analysis.elev{nr}(:,framecount,n) = elev(:,n);
            analysis.energy{nr}(:,framecount,n) = energy;
            analysis.diff{nr}(:,framecount,n) = diffs(:,n); 
        end 
        
      %  disp('overriden params')
%         diffs(:,:) = 1;
      %  azim(:,:) = 0;
      %  elev(:,:) = 0;
         
        %%% SIRR SYNTHESIS %%% 
        z_00 = zeros(nBins_anl, numSec);
        z_ndiff = zeros(nBins_anl, numSec); 
        if pars.RENDER_DIFFUSE
            z_diff = zeros(nBins_anl, numSec); 
        end
        outspec_diff = zeros(nBins_syn, 2);
        outspec_ndiff = zeros(nBins_syn, 2);


        W_S = pars.sectorCoeffs;
        for n=1:numSec  
             
            % NON-DIFFUSE PART
            %diffs_interp = pchip(pars.centerfreqs_anl, diffs(:,n), pars.hrtf_centerfreqs);
            %ndiffs_sqrt = sqrt(1-diffs_interp);
            
            % TODO: handle sqrt(10e-14)
            %ndiffs_sqrt(ndiffs_sqrt<10e-7) = 0;
            
            % Gain factor computation
            %eleindex = round((elev(:,n)*180/pi+90)/vbap_gtable_res(2));
            %aziindex = round(mod(azim(:,n)*180/pi+180,360)/vbap_gtable_res(1));
            %index = aziindex + (eleindex*181) + 1;
            %gains = gtable(index,:);  
            
            % Interpolate azi and ele
            [x, y, z] = sph2cart(azim(:,n), elev(:,n), 1);
            x_inp = pchip(pars.centerfreqs_anl, x, pars.hrtf_centerfreqs);
            y_inp = pchip(pars.centerfreqs_anl, y, pars.hrtf_centerfreqs);
            z_inp = pchip(pars.centerfreqs_anl, z, pars.hrtf_centerfreqs);
            [azim_inp, elev_inp, r_inp] = cart2sph(x_inp, y_inp, z_inp);

            %hrtf_interp = interpHRTFs(rad2deg(azim_inp), rad2deg(elev_inp), pars);
            hrtf_ndif = nearestHRTFs(rad2deg(azim_inp), rad2deg(elev_inp), pars);
            
            %ndiffs_sqrt_inp = pchip(pars.centerfreqs_anl, ndiffs_sqrt, pars.hrtf_centerfreqs);
            % apply ndiff gain to hrtf
            if pars.RENDER_DIFFUSE
                ndiffgains = sqrt(1-diffs(:,n));
                %ndiffgains = (1-diffs(:,n));
            else
                ndiffgains = 1;
            end
            
            % Interpolate panning filters 
            %ndiffgains = interpolateFilters(permute(ndiffgains, [3 2 1]), fftsize);
            %ndiffgains = permute(ndiffgains, [3 2 1]);
            
            % generate non-diffuse stream
            z_00(:,n) = inspec*W_S(:, 4*(n-1) + 1);

            z_ndiff(:, n) = ndiffgains .* z_00(:,n); 
            
            % prepare for frequency domain convolution
            %hrtf_interp = squeeze(permute(interpolateFilters(permute(hrtf_ndif, [3, 2, 1]), fftsize_syn), [3, 2, 1]));
            
            % apply filter
            outspec_ndiff_sec_l = beta_A * interpolateSpectrum(z_ndiff(:, n),fftsize_syn) .* interpolateSpectrum(hrtf_ndif(:,1), fftsize_syn);
            outspec_ndiff_sec_r = beta_A * interpolateSpectrum(z_ndiff(:, n),fftsize_syn) .* interpolateSpectrum(hrtf_ndif(:,2), fftsize_syn);

            outspec_ndiff(:,1) = outspec_ndiff(:,1) + outspec_ndiff_sec_l;
            outspec_ndiff(:,2) = outspec_ndiff(:,2) + outspec_ndiff_sec_r;
            % DIFFUSE PART
            switch pars.RENDER_DIFFUSE
                case 0
                    % No diffuse-field rendering
                case 1
                    % New SIRR diffuse stream rendering, based on re-encoding the 
                    % sector signals scaled with the diffuseness estimates
                    diffgains = sqrt(diffs(:,n));  % TODO use interp
                    %diffgains = (diffs(:,n));  % TODO use interp

                    %diffgains = interpolateFilters(permute(diffgains, [3 2 1]), fftsize);
                    %diffgains = permute(diffgains, [3 2 1]); 
                    if pars.order > 1
                        z_diff(:, n) = diffgains .* z_00(:,n); 
                    end
            end  
         
        if pars.RENDER_DIFFUSE
            %if pars.order == 1
            %    a_diff = repmat(diffgains, [1 nSH]).*inspec;  % Check
            %else
            %    a_diff = z_diff * (beta_A .* Y_enc).';  % Perfect Reconstuction
            %end % encode 
            % prepare for frequency domain convolution
            %a_diff = interpolateSpectrum(a_diff, fftsize_syn);

            %for k=1:nBins_syn
            %    outspec_diff(k,:) = (squeeze(D_bin_interp(:,:,k)) * a_diff(k,:).').'; % decode
            %end
            outspec_diff_sec_l = beta_E * interpolateSpectrum(z_diff(:, n),fftsize_syn) .* interpolateSpectrum(hrtfs_sec(:,1,n), fftsize_syn);
            outspec_diff_sec_r = beta_E * interpolateSpectrum(z_diff(:, n),fftsize_syn) .* interpolateSpectrum(hrtfs_sec(:,2,n), fftsize_syn);

            outspec_diff(:,1) = outspec_diff(:,1) + outspec_diff_sec_l;
            outspec_diff(:,2) = outspec_diff(:,2) + outspec_diff_sec_r;

        end


        end  % sec
        
        % decorrelation based on randomising the phase
        if isequal(pars.decorrelationType, 'phase')
            randomPhi = rand(size(outspec_diff))*2*pi-pi;
            outspec_diff = abs(outspec_diff) .* exp(1i*randomPhi);
        end  
        
        analysis.sf_energy{nr}(framecount,1) = 1/numel(inspec)*trace(inspec'*inspec); 
        analysis.ndiff_energy{nr}(framecount,1) = 1/numel(outspec_ndiff)*trace(outspec_ndiff'*outspec_ndiff); 
        analysis.diff_energy{nr}(framecount,1) = 1/numel(outspec_diff)*trace(outspec_diff'*outspec_diff);
         
        % overlap-add synthesis
        lsir_win_ndiff = real(ifft([outspec_ndiff; conj(outspec_ndiff(end-1:-1:2,:))]));
        lsir_res_ndiff(idx+(0:fftsize_syn-1),:,nr) = lsir_res_ndiff(idx+(0:fftsize_syn-1),:,nr) + lsir_win_ndiff;
        if pars.RENDER_DIFFUSE ~= 0
            lsir_win_diff = real(ifft([outspec_diff; conj(outspec_diff(end-1:-1:2,:))]));
            lsir_res_diff(idx+(0:fftsize_syn-1),:,nr) = lsir_res_diff(idx+(0:fftsize_syn-1),:,nr) + lsir_win_diff;
        end
        
        % advance sample pointer
        idx = idx + hopsize;
        framecount = framecount + 1;
        if framecount >= floor(nFrames/10*progress)  
            fprintf('*');
            progress=progress+1; 
        end  
    end  % //while
    fprintf('\ndone\n')
    
    % remove delay caused by the filter interpolation of gains and circular shift
    delay = fftsize_syn/4+1;  % TODO
    tempout = zeros(size(lsir_res_ndiff(:,:,nr)));
    tempout(1:end-delay,:) = lsir_res_ndiff(1+delay:end,:,nr);
    lsir_res_ndiff(:,:,nr) = tempout;
    if pars.RENDER_DIFFUSE
        tempout = zeros(size(lsir_res_diff(:,:,nr)));
        tempout(1:end-delay,:) = lsir_res_diff(1+delay:end,:,nr);
        lsir_res_diff(:,:,nr) = tempout;
    end
end

% Sum signals at different frequency regions
lsir_pad_ndiff = sum(lsir_res_ndiff, 3);

% Remove delay caused by processing
delay = maxWinsize/2 + xover_order/2; % remove also delay of band-pass filtering
lsir_ndiff = lsir_pad_ndiff(delay + (1:lSig), :);
if pars.RENDER_DIFFUSE
    lsir_pad_diff = sum(lsir_res_diff, 3);
    % Remove delay caused by processing
    lsir_diff = lsir_pad_diff(delay + (1:lSig), :);
end

% apply convolution decorrelation to diffuse stream if specified
if isequal(pars.decorrelationType, 'noise') && pars.RENDER_DIFFUSE
    % we want to apply just enough noise-based reverberation as 
    % to suitably decorrelate the signals, but not change the captured room 
    % characteristics too much. T60s of a very, very dry room should suffice for
    % this task:
    %t60 = [0.07 0.07 0.06 0.04 0.02 0.01];
    t60 = [0.2 0.2 0.16 0.12 0.09 0.04];
    fc = [125 250 500 1e3 2e3 4e3];
    randmat =  synthesizeNoiseReverb(2, pars.fs, t60, fc, 0);
    % randn is gaussian, normalize (FLATTEN achieves the same)
    rand_en = sum(randmat.^2, 1);
    randmat = randmat ./ sqrt(rand_en);
    % Decorrelation
    lsir_diff = fftfilt(randmat, lsir_diff);
    clear randmat;
end 

if pars.BROADBAND_FIRST_PEAK 
    % re-introduce peak based on broadband analysis 
    shir_direct_WXYZ = [shir_direct(:,1).';
        shir_direct(:,4).'./sqrt(3);
        shir_direct(:,2).'./sqrt(3);
        shir_direct(:,3).'./sqrt(3);].';   
    I = real(repmat(conj(shir_direct_WXYZ(:,1)),[1 3]).*shir_direct_WXYZ(:,2:4));
    I = sum(I);
    [dir_azim, dir_elev] = cart2sph(I(1,1), I(1,2), I(1,3));

    % Gain factor computation
%     eleindex = round((dir_elev*180/pi+90)/vbap_gtable_res(2));
%     aziindex = round(mod(dir_azim*180/pi+180,360)/vbap_gtable_res(1));
%     index = aziindex + (eleindex*181) + 1;
%     dir_gains = gtable(index,:);
%     lsir_ndiff = lsir_ndiff + dir_gains .* (shir_direct(:,1)*ones(1,nLS)); 
    
    % Find nearest hrtf to peak DOA, TODO: interpolate?
    [x_p, y_p, z_p] = sph2cart(dir_azim, dir_elev, 1);
    [x_hrfts, y_hrtfs, z_hrtfs] = sph2cart(pars.hrir_dirs_deg(:, 1)*pi/180,...
                                           pars.hrir_dirs_deg(:, 2)*pi/180, 1);
    doa_proj = [x_hrfts, y_hrtfs, z_hrtfs] * [x_p, y_p, z_p].';
    [~, d_min_k] = max(doa_proj);
    p_hrirs = pars.hrirs(:, :, d_min_k);
    
    p_first = sqrt(4*pi) * shir_direct(:,1);% TODO: scale pressure channel?
    lsir_ndiff = lsir_ndiff + fftfilt(p_hrirs, p_first);  % lots of zeros at the end, no padding
end 

if pars.RENDER_DIFFUSE
    lsir = lsir_ndiff + lsir_diff;
else
    lsir = lsir_ndiff;
    lsir_diff = 0;
end 
 
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function M_interp = interpolateFilters(M, fftsize)
% M filter matrix with y = M*x, size NxLxK,
%   N output channels, M input channels, K frequency bins
%
%   Archontis Politis, 12/06/2018
%   archontis.politis@aalto.fi
assert(fftsize >= 2*(size(M, 3)-1))  % CFH: Will truncate otherwise

winsize = 2*(size(M,3)-1);
M_conj = conj(M(:,:,end-1:-1:2));
M_ifft = ifft(cat(3, M, M_conj), [], 3);
M_ifft = M_ifft(:,:, [(winsize/2+1:end) (1:winsize/2)]); % flip
M_interp = fft(M_ifft, fftsize, 3); % interpolate to fftsize
M_interp = M_interp(:,:,1:fftsize/2+1); % keep up to nyquist
end

function X_interp = interpolateSpectrum(X, fftsize)
% Interpolate single sided spectrum X:[t, numCh] by zero padding.
assert(fftsize >= 2*(size(X, 1)-1))  % Will truncate otherwise
assert(mod(size(X, 1), 2) == 1)  % needs to be odd single sided spectrum

X_full = cat(1, X, conj(X(end-1:-1:2, :)));  % mirror
x = ifft(X_full, [], 1);
X_pad = fft(x, fftsize);  % interpolate spectrum by zero padding
X_interp = X_pad(1:fftsize/2+1, :);  % return single sided spectrum
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function sigregs = divide2regions(sig, fs, xovers, fir_order)
% DIVIDE2REGIONS Bandpass the signal for multi-window analysis
%
%   Archontis Politis, 12/06/2018
%   archontis.politis@aalto.fi

% number of frequency regions
nRes = length(xovers)+1;
lIn = size(sig,1);
lOut = lIn+fir_order;
nCH = size(sig,2);
sigregs = zeros(lOut, nCH, nRes);

% create first and last lowpass and highpass in the filterbank
filters = zeros(fir_order+1, nRes);
filters(:,1) = fir1(fir_order, xovers(1)/(fs/2), 'low');
filters(:,nRes) = fir1(fir_order, xovers(nRes-1)/(fs/2), 'high');
for i = 2:(nRes-1)
    filters(:,i) = fir1(fir_order, [xovers(i-1) xovers(i)]/(fs/2), 'bandpass');
end
for i = 1:nRes
    sigregs(:,:,i) = fftfilt(filters(:,i), [sig; zeros(fir_order,nCH)]);
end
%     % omit initial delay of filterbank
%     sigregs = sigregs(fir_order/2+1:end-fir_order/2,:,:);

end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function rir_filt = synthesizeNoiseReverb(nCH, fs, t60, fc, FLATTEN)
%NOISEVERB Simulates a quick and dirty exponential decay reverb tail
%
% order:    HOA order
% fs:       sample rate
% t60:      reverberation times in different bands
% fc:       center frequencies of reverberation time bands (octave bands)
%
%   Archontis Politis, 12/06/2018
%   archontis.politis@aalto.fi

    if nargin<5, FLATTEN = 0; end
    
    % number of HOA channels
    nSH = nCH;
    % number of frequency bands
    nBands = length(t60);
    % decay constants
    alpha = 3*log(10)./t60;
    % length of RIR
    %lFilt = ceil(max(t60)*fs);
    t = (0:1/fs:max(t60)-1/fs)';
    lFilt = length(t);
    % generate envelopes
    env = exp(-t*alpha);
    % generate RIRs
    rir = randn(lFilt, nSH, nBands);
    for k = 1:nBands
        rir(:, :, k) = rir(:,:,k).*(env(:,k)*ones(1,nSH));
    end
    % get filterbank IRs for each band
    filterOrder = 200;
    h_filt = filterbank(fc, filterOrder, fs);
    % filter rirs
    rir_filt = zeros(lFilt+ceil(filterOrder/2), nSH);
    for n = 1:nSH
        h_temp = [squeeze(rir(:,n,:)); zeros(ceil(filterOrder/2), nBands)];
        rir_filt(:, n) = sum(fftfilt(h_filt, h_temp), 2);
    end
                            
    if FLATTEN, rir_filt = equalizeMinphase(rir_filt); end
    
    rir_filt = rir_filt(filterOrder/2+1:end,:); % remove delay
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function h_filt = filterbank(fc, filterOrder, fs)
% fc:   the center frequencies of the bands
% Nord: order of hte FIR filter
%
%   Archontis Politis, 12/06/2018
%   archontis.politis@aalto.fi

    if length(fc) == 1
        h_filt = 1;

    elseif length(fc) == 2
        h_filt = zeros(filterOrder+1, 2);

        % lowpass
        f_ll = 2*fc(1)/sqrt(2);
        w_ll = f_ll/(fs/2);
        h_filt(:, 1) = fir1(filterOrder, w_ll);
        % highpass
        f_hh = fc(2)/sqrt(2);
        w_hh = f_hh/(fs/2);
        h_filt(:, 2) = fir1(filterOrder, w_hh, 'high');

    else
        Nbands = length(fc);
        h_filt = zeros(filterOrder+1, Nbands);

        % lowpass
        f_ll = 2*fc(1)/sqrt(2);
        w_ll = f_ll/(fs/2);
        h_filt(:, 1) = fir1(filterOrder, w_ll);
        % highpass
        f_hh = fc(end)/sqrt(2);
        w_hh = f_hh/(fs/2);
        h_filt(:, end) = fir1(filterOrder, w_hh, 'high');
        % bandpass
        for k = 2:Nbands-1
            fl = fc(k)/sqrt(2);
            fh = 2*fc(k)/sqrt(2);
            wl = fl/(fs/2);
            wh = fh/(fs/2);
            w = [wl wh];
            h_filt(:, k) = fir1(filterOrder, w, 'bandpass');
        end
    end

end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function rir_filt_flat = equalizeMinphase(rir_filt)
%MAKEFLATVERB Makes the decaying noise spectrally flat
%
%   Archontis Politis, 12/06/2018
%   archontis.politis@aalto.fi

Nrir = size(rir_filt,2);
for n=1:Nrir
    % equalise TDI by its minimum phase form to unity magnitude response
    tdi_f = fft(rir_filt(:,n));
    tdi_min_f = exp(conj(hilbert(log(abs(tdi_f)))));
    tdi_eq = real(ifft(tdi_f./tdi_min_f));
    rir_filt_flat(:,n) = tdi_eq;
end

end


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function [sectorCoeffs, secNorms, sec_dirs] = computeSectorCoeffs(orderSec, A_xyz, pattern, grid, sec_dirs)

orderVel = orderSec+1;

if orderSec == 0
    % for N=0, do 1 sector, basic DirAC
    wxyzCoeffs = [ ...
    sqrt(4*pi)      0           0               0;
         0          0           sqrt(4*pi/3)    0;
         0          0           0               sqrt(4*pi/3);
         0       sqrt(4*pi/3)   0               0];

    % convert to real SH coefficients to use with the real signals
    secNorms = [1 1];
    sectorCoeffs = wxyzCoeffs;
    sec_dirs = [];  % or [0, 0], empty more explicit

else
    
    wxyzCoeffs = [];
    if nargin<5
        switch grid
            case 'AP'
                [~, sec_dirs] = getTdesign(orderSec+1); 
            case 'EP'
                [~, sec_dirs] = getTdesign(2*orderSec);
            otherwise
                error("Use 'AP' or 'EP' grid");
        end   
    end
    numSec = size(sec_dirs,1);
    
    switch pattern
        case 'cardioid'
            b_n = beamWeightsCardioid2Spherical(orderSec);
        case 'maxRE'
            b_n = beamWeightsMaxEV(orderSec);
        case {'pwd', 'hypercardioid'}
            b_n = beamWeightsHypercardioid2Spherical(orderSec);
        otherwise
            error("Unknown sector design! " + pattern);
    end
    
    c_n = b_n ./ sqrt((2*(0:orderSec).'+1) / (4*pi));  % remove m
    w_nm = diag(replicatePerOrder(c_n)) * getRSH(orderSec, [0,0]);
    % amplitude preservation factor
    beta_A = (sqrt(4*pi)) / (w_nm(1) * numSec);
    % energy preservation factor (sqrt when applied to signals)
    beta_E = (4*pi) / (w_nm' * w_nm * numSec);

    secNorms = [beta_A, beta_E];
    
    for ns = 1:numSec
        
        % rotate the pattern by rotating the coefficients
        azi_sec = sec_dirs(ns, 1);
        polar_sec = pi/2-sec_dirs(ns, 2); % from elevation to inclination
        c_nm = rotateAxisCoeffs(b_n, polar_sec, azi_sec, 'complex');  % CFH: no compensation here
        % get the velocity coeffs
        x_nm = A_xyz(1:(orderVel+1)^2, 1:(orderSec+1)^2, 1)*c_nm;
        y_nm = A_xyz(1:(orderVel+1)^2, 1:(orderSec+1)^2, 2)*c_nm;
        z_nm = A_xyz(1:(orderVel+1)^2, 1:(orderSec+1)^2, 3)*c_nm;
        
        % Pad the (lower order) sector coefficients and stack in
        % a matrix together with the velocity ones
        c_nm = [c_nm; zeros(2*(orderSec+1)+1,1)];
        wxyzCoeffs = [wxyzCoeffs c_nm x_nm y_nm z_nm];
    end
    % convert to real SH coefficients to use with the real signals
    sectorCoeffs = real(complex2realCoeffs(wxyzCoeffs));
end


end
 
function hrtf_interp = interpHRTFs(azi, elev, pars, freq_bins)
% Azi, Ele in deg
    if nargin < 4
        freq_bins = pars.hrtf_centerfreqs;
    end
    nBins = length(freq_bins);
    assert(length(azi) == nBins);
    hrtf_interp = zeros(nBins, 2);
    
    cutoff = 1500;
    [~,kk_cutoff] = min(abs(freq_bins-cutoff));

    % find closest pre-computed VBAP direction
    az_res = pars.hrtf_vbapTableRes(1); el_res = pars.hrtf_vbapTableRes(2);
    N_azi = round(360/az_res) + 1;
    aziIndex = round(mod(azi+180,360)/double(az_res));
    elevIndex = round((elev+90)/double(el_res));
    gridIndex = elevIndex*double(N_azi)+aziIndex+1; % + 1 for matlab only
    idx3 = pars.hrtf_vbapTableIdx(gridIndex,:);  
    weights3 = pars.hrtf_vbapTableComp(gridIndex,:);
 
    for k=1:nBins
        mags3_nd = reshape(pars.hrtf_mag(k,:,idx3(k,:)),[2 3]).';
        itds3_nd = pars.hrtf_itd(idx3(k,:));
        itd_interp = weights3(k,:) * itds3_nd; 
        mags_interp = weights3(k,:) * mags3_nd;
        
        if k < kk_cutoff
            % convert ITDs to phase differences -pi~pi
            ipd_interp = mod(2*pi*freq_bins(k)*itd_interp + pi, 2*pi) - pi;
        else
            ipd_interp = 0;
        end
        
        % introducing IPDs
        hrtf_interp(k,1) = mags_interp(1,1) .* exp(1i*ipd_interp/2);
        hrtf_interp(k,2) = mags_interp(1,2) .* exp(-1i*ipd_interp/2);
    end
end


function hrtf_nearest = nearestHRTFs(azi, ele, pars, freq_bins)
% Azi, Ele in deg
    if nargin < 4
        freq_bins = pars.hrtf_centerfreqs;
    end
    nBins = length(freq_bins);
    assert(length(azi) == nBins);
    hrtf_nearest = zeros(nBins, 2);

    cutoff = 1500;
    [~,kk_cutoff] = min(abs(freq_bins-cutoff));

    [x_ls, y_ls, z_ls] = sph2cart(azi*pi/180,...
                                  ele*pi/180, 1);
    [x_hrfts, y_hrtfs, z_hrtfs] = sph2cart(pars.hrir_dirs_deg(:, 1)*pi/180,...
                                           pars.hrir_dirs_deg(:, 2)*pi/180, 1);
    ls_proj = [x_hrfts, y_hrtfs, z_hrtfs] * [x_ls, y_ls, z_ls].';
    [~, d_min_k] = max(ls_proj);
%     for k = 1:nBins
%        hrtf_nearest(k, :) = pars.hrtfs(k, :, d_min_k(k));
%     end
    for k = 1:nBins
        if k < kk_cutoff
            % convert ITDs to phase differences -pi~pi
            ipd = mod(2*pi*freq_bins(k)*pars.hrtf_itd(d_min_k(k)) + pi, 2*pi) - pi;
            hrtf_nearest(k, 1) = exp(1i*ipd/2) .* pars.hrtf_mag(k, 1, d_min_k(k));
            hrtf_nearest(k, 2) = exp(-1i*ipd/2) .* pars.hrtf_mag(k, 2, d_min_k(k));
        else
            hrtf_nearest(k, :) = pars.hrtf_mag(k, :, d_min_k(k));
        end
    end
end




