function y = matrixConvolver(x, H, blocksize)
%MATRIXCONVOLVER Summary of this function goes here
%   Detailed explanation goes here
%
%   Implements  in->Filters_matx->out
%
% H:    lh x nCHin x nCHout matrix of filter IRs
%       the filters define a convolution matrix M(f) of size nCHout x nCHin
%       such that y(f) = M(f)*x(f)

nCHin = size(x,2);
nFiltersIn = size(H,2);
nCHout = size(H,3);
if nCHin ~= nFiltersIn
    error('buha');
end

L = blocksize;
lx = size(x,1); % not needed for real-time
lh = size(H,1);
ly = lh+lx-1; % not needed for real-time

% init
numInBlocks = ceil(lx/L); % not needed for real-time
numOutBlocks = ceil(ly/L); % not needed for real-time
numOvrlpAddBlocks = ceil((L+lh-1)/L);
H_f = fft(H,numOvrlpAddBlocks*L,1); % precompute FFTs of filters
H_f = H_f(1:end/2+1,:,:);

y = zeros(numOutBlocks*L,nCHout); % not needed for real-time
x_pad = zeros(numInBlocks*L,nCHin); % not needed for real-time
x_pad(1:lx,:) = x; % not needed for real-time
ovrlpAddBuffer = zeros(numOvrlpAddBlocks*L,nCHout);
HX_n = zeros(numOvrlpAddBlocks*L, nCHin);

% running loop
for n=1:numInBlocks
    x_n = x_pad((n-1)*L+(1:L),:);
    X_n = fft(x_n,numOvrlpAddBlocks*L);
    X_n = X_n(1:end/2+1,:,:);
    for nch = 1:nCHout
        HX_n(1:end/2+1,:) = H_f(:,:,nch).*X_n;
        HX_n(end/2+2:end,:) = conj(HX_n(end/2:-1:2,:));
        hx_n = ifft(HX_n);
        z_n = sum(hx_n,2);
        ovrlpAddBuffer(1:(numOvrlpAddBlocks-1)*L,nch) = ovrlpAddBuffer(L+(1:(numOvrlpAddBlocks-1)*L),nch);
        ovrlpAddBuffer((numOvrlpAddBlocks-1)*L+(1:L),nch) = zeros(L,1);
        ovrlpAddBuffer(:,nch) = ovrlpAddBuffer(:,nch) + z_n;
    end
    y_n = ovrlpAddBuffer(1:L,:); % this is the multichannel buffer going out!
    y((n-1)*L+(1:L),:) = y_n; % not needed for real-time
end
y(numInBlocks*L+(1:(numOvrlpAddBlocks-1)*L),:) = ovrlpAddBuffer(L+(1:(numOvrlpAddBlocks-1)*L),:); % not needed for real-time (handling the last blocks of convolution, after last input block ends
y = y(1:ly,:); % not needed for real-time

end
