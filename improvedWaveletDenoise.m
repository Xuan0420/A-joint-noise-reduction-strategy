function [denoised_signal] = improvedWaveletDenoise(signal, wavelet, level, threshold, alpha, beta)
% Wavelet thresholding signal denoising function (based on the improved threshold function)
% Input parameters:
%   signal - The signal to be denoised
%   wavelet - The chosen wavelet type
%   level - The number of decomposition levels
%   lambda, alpha, beta - The three parameters of the threshold function
% Output parameters:
%   denoised_signal - The denoised signal

% Wavelet decomposition
[c, l] = wavedec(signal, level, wavelet);

% Calculate the threshold and process
c_t = improved_threshold2(c, threshold, alpha, beta);

% Reconstruct the signal
denoised_signal = waverec(c_t, l, wavelet);

end
