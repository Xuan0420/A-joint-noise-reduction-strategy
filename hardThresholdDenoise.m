function [denoised_signal] = waveletDenoise(signal, wavelet, level, threshold)

% Input parameters:
%   signal - The signal to be denoised
%   wavelet - The type of wavelet used
%   level - The number of decomposition levels
%   threshold - The threshold value for soft thresholding
% Output parameters:
%   denoised_signal - The denoised signal

% Wavelet decomposition
[c, l] = wavedec(signal, level, wavelet);

% Calculate hard thresholding
c_t = wthresh(c, 'h', threshold);

% Reconstruct the signal
denoised_signal = waverec(c_t, l, wavelet);

end