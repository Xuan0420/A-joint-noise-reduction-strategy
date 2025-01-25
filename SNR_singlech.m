function snr = SNR_singlech(I, In)
% Calculate the Signal-to-Noise Ratio (SNR) of a noisy signal
% I is the pure signal
% In is the noisy signal
% The formula for SNR calculation is:
% snr = 10 * log10(Esignal / Enoise)

I = I(:)';                                % Convert the data to a column vector
In = In(:)';
Ps = sum((I).^2);                         % Energy of the signal
Pn = sum((I - In).^2);                    % Energy of the noise
snr = 10 * log10(Ps / Pn);                % Ratio of signal energy to noise energy, then convert to decibels
end
