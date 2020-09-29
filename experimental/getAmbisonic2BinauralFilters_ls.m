function [decMtx, decFilters, C, ild, itd, cll, ambihrtfs, ambihrirs] = getAmbisonic2BinauralFilters_ls(hrtfs, hrtf_dirs_deg, order, fs, weights)
%GETAMBISONIC2BINAURALFILTERS_LS Summary of this function goes here
%   Detailed explanation goes here
% hrtfs 2 x nDirs x nBins

nHRTF = size(hrtfs,2);
nBins = size(hrtfs,3);
if nargin<5, weights = ones(nHRTF,1)/nHRTF; end
if weights == -1
    weights = getVoronoiWeights(hrtf_dirs_deg); end
W = diag(weights);

Y_na = getRSH(order, hrtf_dirs_deg);
% Solution matrix
for kk=1:nBins
    H = hrtfs(:,:,kk);
    B_ls = (Y_na*W*Y_na')\(Y_na*W*H'); 
    decMtx(:,:,kk) = B_ls';
end
% Time-domain filters
if nargout>1
    decFilters = cat(3, decMtx, conj(decMtx(:,:,end-1:-1:2)));
    decFilters = ifft(decFilters,[],3);
    decFilters = permute(decFilters,[3 2 1]);
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
    ambihrirs = permute(ambihrirs, [3 1 2]);
    ambihrirs = real(ifft(ambihrirs,[],1));
    itd = computeITDfromXCorr(ambihrirs, fs);
end

end
