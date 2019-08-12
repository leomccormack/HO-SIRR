function [lsir, lsir_ndiff, lsir_diff, pars, analysis] = HOSIRR(shir, pars, ENABLE_PLOTS)
% Higher-Order Spatial Impulse Response Rendering (HOSIRR)
% --------------------------------------------------------
% Multi-resolution Higher-order Spatial Impulse Response Rendering
% 
% DEPENDENCES
%     Spherical-Harmonic-Transform Matlab library
%         https://github.com/polarch/Spherical-Harmonic-Transform
%     Higher-Order-Ambisonics Matlab library
%         https://github.com/polarch/Higher-Order-Ambisonics
%
% INPUT ARGUMENTS
%     shir                   : spherical harmonic domain impulse response
%                              (ACN/N3D), [signalLength  x (order+1)^2]
%     pars.fs                : sample rate
%     pars.ls_dirs_deg       : LS directions in degrees [azi elev]
%     pars.multires_winsize  : [winsize_low winsize_2 ... winsize_high]
%     pars.multires_xovers   : [xover_low ... xover_high]
%     pars.panningNormCoeff  : {0,1} 0 for normal room, 1 for anechoic
%     pars.RENDER_DIFFUSE    : {0,1} 0 off, 1 new diffuse stream via
%                              ambisonic decoding of diffuseness scaled
%                              sector signals
%     pars.BROADBAND_DIRECT  : {0,1} 0 off, 1 broadband analysis for direct
%     pars.nBroadbandPeaks   : number of peaks to pan using broadband DoAs
%     pars.cycles_diff       : time constant for diffuseness computation, 
%                              in cycles per center frequency
%     pars.alpha_diff        : smallest one-pole allowed alpha value for
%                              smoothing diff 
%                              y(n) = alpha*y(n-1) + (1-alpha)*x(n)
%     pars.decorrelationType : {'phase','noise'}  
%     pars.orderPerBand      : scalar for same order at all bands,
%                              [nBandsx1] for different orders at different
%                              bands (NOT IMPLEMENTED YET)
%
% OUTPUT ARGUMENTS
%     lsir                   : loudspeaker impulse responses    
%     lsir_ndiff             : non-diffuse stream only
%     lsir_diff              : diffuse stream only
%     pars                   : parameter struct used for the processing
%     analysis               : {nRes}[nBins x nFrames x nSectors], historic
%                              DoA estimates (.azim, .elev), .energy and
%                              diffuseness (.diff) estimates
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
%   Leo McCormack, 22/09/2018
%   leo.mccormack@aalto.fi
%   Archontis Politis, 12/06/2018
%   archontis.politis@aalto.fi
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%%% Defaults/Warnings/Errors: 
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

disp('HOSIRR Configuration:'), pars %#ok

lSig = size(shir,1);
nSH = size(shir,2);

% normalise input to max(|insig|) = 1
shir = shir./(max(abs(shir(:,1))));

% extract the first N highest peaks
if pars.BROADBAND_DIRECT
    shir_tmp = shir;
    inv_direct_win = ones(lSig,1);
    figure
    for peak = 1:pars.nBroadbandPeaks
        % find the index of highest peak in the omni
        [~, peak_ind] = max(abs(shir_tmp(:,1)).^2);

        % calculate window 
        dirwinsize = 64;
        direct_win = zeros(lSig,1);
        direct_win(peak_ind-dirwinsize/2:peak_ind+dirwinsize/2,1) = hanning(dirwinsize+1);
        if ENABLE_PLOTS
            plot(abs(shir_tmp(:,1)).^2, 'b'), hold on
            plot(direct_win.*abs(shir_tmp(:,1)).^2, 'r'), hold on
            plot(direct_win, 'r'), title('windowed'), hold on
            xlim([0 lSig/16])
            legend('IR', 'Windowed')
        end 

        % extract peak from shir
        shir_direct{peak} = repmat(direct_win, [1 nSH]).*shir; %#ok   
        shir_tmp = repmat(1-direct_win, [1 nSH]).*shir_tmp; 
        inv_direct_win = (1-direct_win).*inv_direct_win; 
    end
    shir = shir_tmp; % remove peaks for the main loop
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
maxOrder = sqrt(nSH)-1;
for order_i = 1:maxOrder
    [~,sec_dirs_rad] = getTdesign(2*(order_i)); 
    A_xyz = computeVelCoeffsMtx(order_i-1);
    [pars.sectorCoeffs{order_i}, pars.normSec{order_i}] = computeSectorCoeffs(order_i-1, A_xyz, 'pwd', sec_dirs_rad, 'EP');
    if order_i~=1
        pars.sectorDirs{order_i-1} = sec_dirs_rad;
    end
end 
maxNumSec = size(pars.sectorCoeffs{maxOrder},2)/4;

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
lsir_res_ndiff = zeros(lSig + 5*maxfftsize + xover_order, nLS, nRes);
lsir_res_diff = zeros(lSig + 5*maxfftsize + xover_order, nLS, nRes); 

%%% Multi-resolution processing
for nr = 1:nRes 
    disp(['Processing frequency region no. ' num2str(nr)]);
    winsize = pars.multires_winsize(nr);
    fftsize = 2*winsize; % double the window size for FD convolution
    hopsize = winsize/2; % half the window size time-resolution
    nBins_anl = winsize/2 + 1; % nBins used for analysis
    nBins_syn = fftsize/2 + 1; % nBins used for synthesis
    centerfreqs = (0:fftsize/2)'*pars.fs/fftsize;
    if nr==nRes
        [~,analysisLimit_ind] = min(abs(centerfreqs-min(pars.maxDiffuseAnalysis_Hz)));
    else
        [~,analysisLimit_ind] = min(abs(centerfreqs-min(pars.maxDiffuseAnalysis_Hz,pars.multires_xovers(nr))));  
    end 
    
    % storage for estimated parameters
    analysis.azim{nr} = nan(nBins_anl, ceil(lSig/hopsize), maxNumSec);
    analysis.elev{nr} = nan(nBins_anl, ceil(lSig/hopsize), maxNumSec);
    analysis.energy{nr} = nan(nBins_anl, ceil(lSig/hopsize), maxNumSec);
    analysis.diff{nr} = nan(nBins_anl, ceil(lSig/hopsize), maxNumSec);
     
    % extended energy analysis
    analysis.sf_energy{nr} = nan(ceil(lSig/hopsize),1);
    analysis.sec_energy{nr} = nan(ceil(lSig/hopsize),1);
    analysis.ndiff_energy{nr} = nan(ceil(lSig/hopsize),1);
    analysis.diff_energy{nr} = nan(ceil(lSig/hopsize),1);
    
    % transform window (hanning)
    x = 0:(winsize-1);
    win = sin(x.*(pi/winsize))'.^2;
    
    % sectors
    order = maxOrder; 
    sectorCoeffs_order = pars.sectorCoeffs{order}./sqrt(4*pi);  
    nSec_order = size(sectorCoeffs_order,2)/4;
    
    % Loudspeaker to sector assignment (HOA only)
    if order>1
        sectorDirs_order = pars.sectorDirs{order-1}; 
    end 
    
    switch pars.RENDER_DIFFUSE
        case 0
            % No diffuse stream
        case 1
        % New SIRR diffuse stream, based on scaling the sector signals with diffuseness estimates
        % and then re-encoding them into SHs, and then decoding them to the loudspeaker setup
            if order>1
                Y_enc = getRSH(order, sectorDirs_order*180/pi); % encoder
            end   
            M_diff = getRSH(order, ls_dirs_deg).'/sqrt(nLS*nSH);   
    end 
     
    % diffuseness averaging buffers 
    prev_intensity = zeros(3,nSec_order);
    prev_energy = zeros(nSec_order,1);
    
    % analysed parameters for each sector
    azim = zeros(nBins_anl,nSec_order);
    elev = zeros(nBins_anl,nSec_order);
    diffs = zeros(nBins_anl,nSec_order); 
    
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
        WXYZ_ana = inspec_anl(:,1:(order+1)^2)*sectorCoeffs_order; 
        for n=1:nSec_order 
            % form weighted pressure-velocity signals
            WXYZ_sec = WXYZ_ana(:,4*(n-1) + (1:4)); 
            energy = 0.5.*sum(abs(WXYZ_sec).^2,2);
            
            % Frequency-dependent instantanious DoA estimation
            I = real(conj(WXYZ_sec(:,1)*ones(1,3)) .* WXYZ_sec(:,2:4));
            [azim(:,n), elev(:,n)] = cart2sph(I(:,1), I(:,2), I(:,3));
            
            % Broad-band active intensity vector for diffuseness estimation
            pvCOV = (WXYZ_sec(1:analysisLimit_ind,:)'*WXYZ_sec(1:analysisLimit_ind,:)); 
            Ia = real(pvCOV(2:4,1)); 
            E = trace(pvCOV)/2;  
            
            % Time averaging of boadband diffuseness
            diff_intensity = (1-pars.alpha_diff).*Ia + pars.alpha_diff.*prev_intensity(:,n);
            diff_energy = (1-pars.alpha_diff).*E + pars.alpha_diff.*prev_energy(n); 
            prev_intensity(:,n) = diff_intensity;
            prev_energy(n) = diff_energy; 
            diffs(:,n) = 1 - sqrt(sum(diff_intensity.^2)) ./ (diff_energy + eps); 
             
            % storage for estimated parameters over time
            analysis.azim{nr}(:,framecount,n) = azim(:,n);
            analysis.elev{nr}(:,framecount,n) = elev(:,n);
            analysis.energy{nr}(:,framecount,n) = energy;
            analysis.diff{nr}(:,framecount,n) = diffs(:,n); 
        end 
        
        %%% SIRR SYNTHESIS %%% 
        if pars.RENDER_DIFFUSE, W_diff = zeros(nBins_syn, nSec_order); end
        W_syn = zeros(nBins_syn, nSec_order);
        for n=1:nSec_order  
             
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

            % Normalisation term
            nnorm = sqrt(pars.normSec{order});
            
            % generate non-diffuse stream
            W_syn(:,n) = inspec_syn(:,1:(order+1)^2)*sectorCoeffs_order(:, 4*(n-1) + 1);
            outspec_ndiff = outspec_ndiff + ndiffgains .* (nnorm.*W_syn(:,n)*ones(1,nLS));
    
            % DIFFUSE PART
            switch pars.RENDER_DIFFUSE
                case 0
                    % No diffuse-field rendering
                case 1
                    % New SIRR diffuse stream rendering, based on re-encoding the 
                    % sector signals scaled with the diffuseness estimates
                    diffgains = sqrt(diffs(:,n));  
                    diffgains = interpolateFilters(permute(diffgains, [3 2 1]), fftsize);
                    diffgains = permute(diffgains, [3 2 1]); 
                    if order == 1, a_diff = repmat(diffgains, [1 nSH]).*inspec_syn./sqrt(nSH);
                    else, W_diff(:, n) = diffgains .* W_syn(:,n); end   
            end  
        end 
        if pars.RENDER_DIFFUSE
            if order > 1, a_diff = W_diff./sqrt(nSec_order) * Y_enc.'; end % encode 
            outspec_diff = sqrt(nSH) .* a_diff * M_diff.'; % decode
        end
        
        % decorrelation based on randomising the phase
        if isequal(pars.decorrelationType, 'phase')
            randomPhi = rand(size(outspec_diff))*2*pi-pi;
            outspec_diff = abs(outspec_diff) .* exp(1i*randomPhi);
        end  
        
        analysis.sf_energy{nr}(framecount,1) = mean(sum(abs(inspec_syn).^2/nSH,2));
        analysis.sec_energy{nr}(framecount,1) = mean(sum(abs(W_syn).^2,2));
        analysis.ndiff_energy{nr}(framecount,1) = mean(sum(abs(outspec_ndiff).^2,2));
        analysis.diff_energy{nr}(framecount,1) = mean(sum(abs(outspec_diff).^2,2));
         
%         % ambi_ = mean(sum(abs(a_diff).^2,2)); 
%         sf_ = analysis.sf_energy{nr}(framecount,1);
%         sec_ = analysis.sec_energy{nr}(framecount,1);
%         ndiff_ = analysis.ndiff_energy{nr}(framecount,1);
%         diff_ = analysis.diff_energy{nr}(framecount,1);
%          
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
    % Decorrelation matrix 
    t60 = [0.07 0.07 0.06 0.04 0.02 0.01];
    fc = [125 250 500 1e3 2e3 4e3];
    randmat =  synthesizeNoiseReverb(nLS, pars.fs, t60, fc, 1); 
    % Decorrelation
    lsir_diff = fftfilt(randmat, lsir_diff);
    clear randmat;
end 

if pars.BROADBAND_DIRECT
    % remove peaks from non-diffuse stream
    %lsir_ndiff = lsir_ndiff.*repmat(inv_direct_win, [1, nLS]); 
    %ls_energy = zeros(lSig, 1);
    
    % re-introduce peaks based on broadband analysis
    for peak = 1:pars.nBroadbandPeaks 
        shir_direct_WXYZ = [shir_direct{peak}(:,1).';
            shir_direct{peak}(:,4).'./sqrt(3);
            shir_direct{peak}(:,2).'./sqrt(3);
            shir_direct{peak}(:,3).'./sqrt(3);].';   
        I = real(repmat(conj(shir_direct_WXYZ(:,1)),[1 3]).*shir_direct_WXYZ(:,2:4));
        I = sum(I); %sum(I.* repmat(abs(shir_direct_WXYZ(:,1)), [1 3])); % weighted average
        [dir_azim, dir_elev] = cart2sph(I(1,1), I(1,2), I(1,3));

        % Gain factor computation
        eleindex = round((dir_elev*180/pi+90)/vbap_gtable_res(2));
        aziindex = round(mod(dir_azim*180/pi+180,360)/vbap_gtable_res(1));
        index = aziindex + (eleindex*181) + 1;
        dir_gains = gtable(index,:);

        % Add the first peak
        %ls_energy = ls_energy + sum( abs(dir_gains .* (shir_direct{peak}(:,1)*ones(1,nLS))).^2, 2);
        lsir_ndiff = lsir_ndiff + dir_gains .* (shir_direct{peak}(:,1)*ones(1,nLS));
    end 
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
function [sectorCoeffs, normSec] = computeSectorCoeffs(orderSec, A_xyz, pattern, sec_dirs, norm)
% COMPUTESECTORCOEFFS Computes the beamforming matrices of sector and 
% velocity coefficients for energy-preserving sectors, for orderSec, 
% for real SH.
%
%   Archontis Politis, 12/06/2018
%   archontis.politis@aalto.fi

orderVel = orderSec+1;

if orderSec == 0
    % for N=0, do 1 sector, basic DirAC
    wxyzCoeffs = [ ...
    sqrt(4*pi)      0           0               0;
         0          0           sqrt(4*pi/3)    0;
         0          0           0               sqrt(4*pi/3);
         0       sqrt(4*pi/3)   0               0];

    % convert to real SH coefficients to use with the real signals
    normSec = 1;
    sectorCoeffs = normSec*wxyzCoeffs;
    
else
    
    wxyzCoeffs = [];
    if nargin<4
        switch norm
            case 'AP'
                [~, sec_dirs] = getTdesign(orderSec+1); 
            case 'EP'
                [~, sec_dirs] = getTdesign(2*orderSec);
        end 
        
    end
    numSec = size(sec_dirs,1);
    
    switch pattern
        case 'cardioid'
            b_n = beamWeightsCardioid2Spherical(orderSec);
            Q = 2*orderSec+1;
        case 'maxRE'
            b_n = beamWeightsMaxEV(orderSec);
            Q = 4*pi/(b_n'*b_n);            
        case 'pwd'
            b_n = beamWeightsHypercardioid2Spherical(orderSec);
            Q = (orderSec+1)^2;            
    end
    
    switch norm
        case 'AP'
            % amplitude normalisation for sector patterns
            normSec = (orderSec+1)/numSec;
        case 'EP'
            % energy normalisation for sector patterns
            normSec = Q/numSec;
    end 
    
    for ns = 1:numSec
        
        % rotate the pattern by rotating the coefficients
        azi_sec = sec_dirs(ns, 1);
        polar_sec = pi/2-sec_dirs(ns, 2); % from elevation to inclination
        c_nm = sqrt(normSec) * rotateAxisCoeffs(b_n, polar_sec, azi_sec, 'complex');
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
 








