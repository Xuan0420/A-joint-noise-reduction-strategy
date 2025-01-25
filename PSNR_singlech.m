function psnr = PSNR_singlech(I, In)
% Calculate the Peak Signal-to-Noise Ratio (PSNR) of a noisy signal
% I is the pure signal
% In is the noisy signal
% The formula for PSNR calculation is:
% psnr = 10 * log10(peak^2 / MSE)

I = I(:)';                                  % Convert the data to a column vector
In = In(:)';
MSE = mean((I - In).^2);                    % Calculate the Mean Squared Error (MSE)
peak = max(I);                              % Find the peak value (maximum amplitude)
psnr = 10 * log10(peak.^2 / MSE);            % Calculate the Peak Signal-to-Noise Ratio in decibels
end
