function [snr, psnr, rmse] = EvaMetrix(x, xwpd)

% Denoising effect is evaluated based on Mean Squared Error and Signal-to-Noise Ratio
snr = SNR_singlech(x, xwpd);        % Calculate Signal-to-Noise Ratio
psnr = PSNR_singlech(x, xwpd);      % Calculate Peak Signal-to-Noise Ratio
rmse = sqrt(sum((xwpd - x).^2)) / length(x); % Calculate Root Mean Squared Error
% Calculate signal distortion
% disp(['Filtered Signal-to-Noise Ratio SNR: ' num2str(snr)])
% disp(['Filtered Root Mean Squared Error RMSE: ' num2str(rmse)])
