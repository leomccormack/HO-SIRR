function [decMtx, decFilters, C, ild, itd, cll, ambihrtfs, ambihrirs] = getAmbisonic2BinauralFilters_magls_zotter(hrtfs, hrtf_dirs_deg, order, fs, cutoff, weights)
%GETAMBISONIC2BINAURALFILTERS_LS Summary of this function goes here
%   Detailed explanation goes here
% hrtfs 2 x nDirs x nBins
% decMtx 2 x nSH x nBins

if nargin<6, weights = false; end  % should sum to 4pi
if nargin<5, cutoff = 1500; end

nHRTF = size(hrtfs,2);
nBins = size(hrtfs,3);
nSH = (order+1)^2;

W = diag(weights);

Y_na = getRSH(order, hrtf_dirs_deg);

% Remove itd from high frequency HRTFs
nFFT = 2*(nBins-1);
f = (0:nFFT/2)*fs/nFFT;
[~,kk_cutoff] = min(abs(f-cutoff));
% Solution matrix
B_magls = zeros(nSH, 2, nBins);
for kk=1:nBins
    if kk<kk_cutoff 
        H = hrtfs(:,:,kk);
        if weights
            B_magls(:,:,kk) = (Y_na*W*Y_na.')\(Y_na*W*H.');
        else
            %B_magls(:,:,kk) = pinv(Y_na.') * H.';
            B_magls(:,:,kk) = (Y_na*Y_na.')\Y_na*H.';
        end
    else
%         D_prev = decMtx(:,:,kk-1);%B_magls(:,:,kk-1)';
%         H_ambi_prev = D_prev*Y_na;
%         phiH_prev = angle(H_ambi_prev);
%         H_mod = abs(hrtfs(:,:,kk)).*exp(1i*phiH_prev);
%         B_magls(:,:,kk) = (Y_na*W*Y_na')\(Y_na*W*H_mod.');
        
        D_prev = B_magls(:,:,kk-1);
        H_ambi_prev = Y_na.' * D_prev;
        phiH_prev = angle(H_ambi_prev);
        H_mod = abs(hrtfs(:,:,kk)).*exp(1i*phiH_prev.');
        if weights
            B_magls(:,:,kk) = (Y_na*W*Y_na.')\(Y_na*W*H_mod.');
        else
            %B_magls(:,:,kk) = pinv(Y_na.') * H_mod.';
            B_magls(:,:,kk) = (Y_na*Y_na.')\Y_na*H_mod.';
        end
    end
end
%decMtx = conj(permute(B_magls, [2, 1, 3]));  % B_magls(:,:,kk)'
decMtx = permute(B_magls, [2, 1, 3]);

% Time-domain filters
if nargout>1
    decFilters = cat(3, decMtx, conj(decMtx(:,:,end-1:-1:2)));
    decFilters = real(ifft(decFilters,[],3));
    %decFilters = permute(decFilters,[3 2 1]);
end
% Reconstructed hrtfs, diffuse-field coherence matrix, ITDs, ILDs
if nargout>2
    for kk=1:nBins
        ambihrtfs(:,:,kk) = decMtx(:,:,kk)*Y_na;
        C(:,:,kk) = ambihrtfs(:,:,kk)*W*ambihrtfs(:,:,kk)';
        ild(:,kk) = 10*log10(abs(ambihrtfs(1,:,kk)).^2  ./ abs(ambihrtfs(2,:,kk)).^2);
        cll(:,kk) = 10*log10(sum( abs(ambihrtfs(:,:,kk)).^2 ) ./ sum( abs(hrtfs(:,:,kk)).^2 ) );
    end
    ambihrirs = cat(3, ambihrtfs, conj(ambihrtfs(:,:,end-1:-1:2)));
    ambihrirs = real(ifft(ambihrirs,[],3));
    itd = computeITDfromXCorr(permute(ambihrirs, [3 1 2]), fs);
end

end
