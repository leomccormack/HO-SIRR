function itds = computeITDfromXCorr(hrirs, fs, cutoff)
%COMPUTEITDFROMXCORR Summary of this function goes here
%   Detailed explanation goes here

if nargin<3, cutoff = 800; end

nDirs = size(hrirs,3);

lp = fir1(250, cutoff/(fs/2),'low');
if cutoff
    hrirs_lp = filter(lp,1,[hrirs; zeros(250,2,nDirs)]);
else
    hrirs_lp = hrirs;
end

lHrir = size(resample(hrirs_lp(:,1,1),8*fs,fs),1);
% itds measured from as delay from left ear first - positive delay means
% source on the left hemisphere (right hrir after the left), negative delay 
% means source on the right (right hrir before left)
for n=1:nDirs
    h = squeeze(hrirs_lp(:,:,n));
    h = resample(h,8*fs,fs);
    [~, maxcorr] = max(xcorr(h(:,1), h(:,2)));
    itds(n,1) = (lHrir - maxcorr)/(8*fs);
end

end
