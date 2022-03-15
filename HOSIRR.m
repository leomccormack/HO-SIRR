function [lsir, lsir_ndiff, lsir_diff, pars, analysis] = HOSIRR(shir, pars)
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
%   Leo McCormack, 22/09/2018
%   leo.mccormack@aalto.fi
%   Archontis Politis, 12/06/2018
%   archontis.politis@aalto.fi
%   Chris Hold, 17/11/2021
%   christoph.hold@aalto.fi
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
if (~isfield(pars, 'ls_dirs_deg')) 
    error('Please specify "ls_dirs_deg", in degrees')
end 
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

% Loudspeaker locations and VBAP gain table
ls_dirs_deg = pars.ls_dirs_deg;
nLS = size(ls_dirs_deg,1);
vbap_gtable_res = [2 2]; % azi, elev, step sizes in degrees
gtable = getGainTable(ls_dirs_deg, vbap_gtable_res); 
 
% Sector design
if (~isfield(pars, 'pattern')), pars.pattern = 'pwd'; end
if (~isfield(pars, 'ENABLE_DIFF_PR')), pars.ENABLE_DIFF_PR = true; end

[~,sec_dirs_rad] = getTdesign(2*(pars.order));
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
maxfftsize = 2*maxWinsize; 
lsir_res_ndiff = zeros(lSig + 2*maxfftsize + xover_order, nLS, nRes);
lsir_res_diff = zeros(lSig + 2*maxfftsize + xover_order, nLS, nRes); 

%%% Multi-resolution processing
for nr = 1:nRes 
    disp(['Processing frequency region no. ' num2str(nr)]);
    winsize = pars.multires_winsize(nr);
    fftsize = 2*winsize; % double the window size for FD convolution
    hopsize = winsize/2; % half the window size time-resolution
    nBins_anl = winsize/2 + 1; % nBins used for analysis
    nBins_syn = fftsize/2 + 1; % nBins used for synthesis 
    centerfreqs_anl = (0:winsize/2)'*pars.fs/winsize;
    if pars.BROADBAND_DIFFUSENESS
        if nr==nRes
            [~,maxDiffFreq_Ind] = min(abs(centerfreqs_anl-min(pars.maxDiffFreq_Hz)));
        else
            [~,maxDiffFreq_Ind] = min(abs(centerfreqs_anl-min(pars.maxDiffFreq_Hz,pars.multires_xovers(nr))));  
        end 
    end  
    
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
            % SHs, and then decoding them to the loudspeaker setup
            if pars.order==1
                D_ls = sqrt(4*pi/nLS).*getRSH(pars.order, ls_dirs_deg).';   
            else 
                Y_enc = diag(replicatePerOrder(c_n_res))*...
                    getRSH(pars.order-1, pars.sectorDirs*180/pi); % encoder
                D_ls = sqrt(4*pi/nLS).*getRSH(pars.order-1, ls_dirs_deg).';   
            end    
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
        inspec_syn = fft(insig_win, fftsize);
        inspec_syn = inspec_syn(1:nBins_syn,:); % keep up to nyquist
        
        % Do analysis using only true resolution
        inspec_anl = inspec_syn(1:fftsize/winsize:end,:);
         
        %%% SIRR ANALYSIS %%%
        outspec_ndiff = 0;
        outspec_diff = 0;
        W_S = pars.sectorCoeffs;   
        s_ana = inspec_anl*W_S; 
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
                energy = 0.5 .* real(trace(pvCOV));  % cast back to real

                % Estimating and time averaging of boadband diffuseness
                diff_intensity = (1-pars.alpha_diff).*I_diff + pars.alpha_diff.*prev_intensity(:,n);
                diff_energy = (1-pars.alpha_diff).*energy + pars.alpha_diff.*prev_energy(n); 
                prev_intensity(:,n) = diff_intensity;
                prev_energy(n) = diff_energy; 
                diffs(:,n) = 1 - sqrt(sum(diff_intensity.^2)) ./ (diff_energy + eps); 
            else  
                energy = 0.5 .* sum(abs(WXYZ_sec).^2,2); 
            
                % Time averaging of intensity-vector for the diffuseness
                % estimate per bin
                diff_intensity = (1-pars.alpha_diff).*I + pars.alpha_diff.*prev_intensity(:,:,n);
                diff_energy = (1-pars.alpha_diff).*energy + pars.alpha_diff.*prev_energy(:,n); 
                diffs(:,n) = 1 - sqrt(sum(diff_intensity.^2,2)) ./ (diff_energy + eps); 
                prev_intensity(:,:,n) = diff_intensity;
                prev_energy(:,n) = diff_energy;
                %assert(all(diffs(:,n)<=1.001))
                %assert(all(diffs(:,n)>=0))
            end 
             
            % storage for estimated parameters over time
            analysis.azim{nr}(:,framecount,n) = azim(:,n);
            analysis.elev{nr}(:,framecount,n) = elev(:,n);
            assert(all(energy>=0))
            analysis.energy{nr}(:,framecount,n) = energy;
            analysis.diff{nr}(:,framecount,n) = diffs(:,n); 
        end 
        
        %%% SIRR SYNTHESIS %%% 
        if pars.RENDER_DIFFUSE
            z_diff = zeros(nBins_syn, numSec); 
        end
        z_00 = zeros(nBins_syn, numSec); 
        W_S = pars.sectorCoeffs;   
        for n=1:numSec  
             
            % NON-DIFFUSE PART
            ndiffs_sqrt = sqrt(1-diffs(:,n)); 
            
            % Gain factor computation
            eleindex = round((elev(:,n)*180/pi+90)/vbap_gtable_res(2));
            aziindex = round(mod(azim(:,n)*180/pi+180,360)/vbap_gtable_res(1));
            index = aziindex + (eleindex*181) + 1;
            gains = gtable(index,:);  
            
            % Interpolate the gains in frequency for proper convolution
            if pars.RENDER_DIFFUSE
                ndiffgains = gains .* (ndiffs_sqrt*ones(1,nLS)); 
            else
                ndiffgains = gains; 
            end
            
            % Interpolate panning filters 
            ndiffgains = interpolateFilters(permute(ndiffgains, [3 2 1]), fftsize);
            ndiffgains = permute(ndiffgains, [3 2 1]);
            
            % generate non-diffuse stream
            z_00(:,n) = inspec_syn*W_S(:, 4*(n-1) + 1);
            outspec_ndiff = outspec_ndiff + ndiffgains .* (beta_A.*z_00(:,n)*ones(1,nLS));
    
            % DIFFUSE PART
            switch pars.RENDER_DIFFUSE
                case 0
                    % No diffuse-fieldendering
                case 1
                    % New SIRR diffuse stream rendering, based on re-encoding the 
                    % sector signals scaled with the diffuseness estimates
                    diffgains = sqrt(diffs(:,n));  
                    diffgains = interpolateFilters(permute(diffgains, [3 2 1]), fftsize);
                    diffgains = permute(diffgains, [3 2 1]); 
                    if pars.order == 1
                        a_diff = repmat(diffgains, [1 nSH]).*inspec_syn;  % Check
                    else
                        z_diff(:, n) = diffgains .* sqrt(beta_E) .* z_00(:,n); 
                    end
            end  
        end 
        if pars.RENDER_DIFFUSE
            if pars.order > 1
                a_diff = z_diff * Y_enc.'; 
            end % encode 
            outspec_diff = a_diff * D_ls.'; % decode
        end
        
        % decorrelation based on randomising the phase
        if isequal(pars.decorrelationType, 'phase')
            randomPhi = rand(size(outspec_diff))*2*pi-pi;
            outspec_diff = abs(outspec_diff) .* exp(1i*randomPhi);
        end  
        
        analysis.sf_energy{nr}(framecount,1) = mean((4*pi/nSH)*sum(abs(inspec_syn).^2,2)); 
        analysis.ndiff_energy{nr}(framecount,1) = mean(sum(abs(outspec_ndiff).^2,2)); 
        analysis.diff_energy{nr}(framecount,1) = mean(sum(abs(outspec_diff).^2,2));
        analysis.total_energy{nr}(framecount,1) = mean(sum(abs(outspec_ndiff+outspec_diff).^2,2));


        % overlap-add synthesis
        lsir_win_ndiff = real(ifft([outspec_ndiff; conj(outspec_ndiff(end-1:-1:2,:))]));
        lsir_res_ndiff(idx+(0:fftsize-1),:,nr) = lsir_res_ndiff(idx+(0:fftsize-1),:,nr) + lsir_win_ndiff;
        if pars.RENDER_DIFFUSE ~= 0
            lsir_win_diff = real(ifft([outspec_diff; conj(outspec_diff(end-1:-1:2,:))]));
            lsir_res_diff(idx+(0:fftsize-1),:,nr) = lsir_res_diff(idx+(0:fftsize-1),:,nr) + lsir_win_diff;
        end
        
        % advance sample pointer
        idx = idx + hopsize;
        framecount = framecount + 1;
        if framecount >= floor(nFrames/10*progress)  
            fprintf('*');
            progress=progress+1; 
        end  
    end
    fprintf('\ndone\n')
    
    % remove delay caused by the filter interpolation of gains and circular shift
    tempout = zeros(size(lsir_res_ndiff(:,:,nr)));
    tempout(1:end-winsize/2,:) = lsir_res_ndiff(winsize/2+1:end,:,nr);
    lsir_res_ndiff(:,:,nr) = tempout;
    if pars.RENDER_DIFFUSE
        tempout = zeros(size(lsir_res_diff(:,:,nr)));
        tempout(1:end-winsize/2,:) = lsir_res_diff(winsize/2+1:end,:,nr);
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
    randmat =  synthesizeNoiseReverb(nLS, pars.fs, t60, fc, 1); 
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
    eleindex = round((dir_elev*180/pi+90)/vbap_gtable_res(2));
    aziindex = round(mod(dir_azim*180/pi+180,360)/vbap_gtable_res(1));
    index = aziindex + (eleindex*181) + 1;
    dir_gains = gtable(index,:);

    lsir_ndiff = lsir_ndiff + dir_gains .* (shir_direct(:,1)*ones(1,nLS)); 
end 

if pars.RENDER_DIFFUSE
    lsir = lsir_ndiff+lsir_diff;
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

winsize = 2*(size(M,3)-1);
M_conj = conj(M(:,:,end-1:-1:2));
M_ifft = ifft(cat(3, M, M_conj), [], 3);
M_ifft = M_ifft(:,:, [(winsize/2+1:end) (1:winsize/2)]); % flip
M_interp = fft(M_ifft, fftsize, 3); % interpolate to fftsize
M_interp = M_interp(:,:,1:fftsize/2+1); % keep up to nyquist
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
 
 
