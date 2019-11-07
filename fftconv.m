function y = fftconv(h, x)
%FFTCONV Linear convolution of x and h through FFT
%
%   FFTCONV(X,H) convolves X and H by multiplying their fourier transforms.
%   Before the FFT transforms the sequences are zero padded to the next 
%   power of 2 higher than length(X)+length(Y)-1 , to avoid circular
%   convolution artifacts. Then an inverse FFT is applied to the product.
%
%	h: IR
% 	x: input signal
%
% 	Written by Archontis Politis, archontis.politis@aalto.fi

lh = length(h);
lx = length(x);
ly = lh+lx-1;

h0 = zeros(2^nextpow2(ly),1);
x0 = zeros(2^nextpow2(ly),1);
h0(1:lh) = h;
x0(1:lx) = x;

y = real(ifft(fft(h0).*fft(x0)));
y = y(1:ly);

end
