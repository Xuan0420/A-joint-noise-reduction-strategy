function rs = EMD(t, fs, y, ys)
%% EMD
u = emd(ys);
[m,~] = size(u);

figure
for i=1:m-1
    subplot(m,1,i)
    plot(t,u(i,:),'b-','linewidth',1)
    ylabel(['IMF',num2str(i)]);
    axis tight
end
subplot(m,1,m)
plot(t,u(m,:),'b-','linewidth',1)
ylabel('Res');
axis tight

figure
for i=1:m-1
    subplot(m,1,i)
    
    %% FFT
    [cc, y_f] = hua_fft(u(i,:), fs, 1);
    plot(y_f, cc, 'b', 'LineWidth', 1.5);
    ylabel(['FFT IMF', num2str(i)]);
    axis tight
end
subplot(m,1,m)
[cc, y_f] = hua_fft(u(m,:), fs, 1);
plot(y_f, cc, 'b', 'LineWidth', 1.5);
axis tight
ylabel('FFT Res');

%% Calculate envelope entropy
[m, n] = size(u);
for i=1:m
    En(i) = BaoluoEntropy(u(i,:));
end
disp(['Envelope entropy of each IMF: ', num2str(En)])

%% HHT spectrogram
[A, f, tt] = hhspectrum(u); % Apply Hilbert transform to each IMF, then apply FFT. Here, f is the normalized frequency.

%% A: instantaneous amplitude, f: instantaneous frequency
%% Synthesize HHT spectrogram
[E, ttt, Cenf] = toimage(A, f);

%% Plot HHT spectrogram
figure
imagesc(t, [0, 0.5*fs], E);
set(gca, 'YDir', 'normal')
text(0.05, 0.9, '(a)', 'Units', 'Normalized', 'FontSize', 14, 'Color', 'white');
xlabel('Time/s')
ylabel('Frequency/Hz')
title('Hilbert spectrum (EMD)')
ylim([0, 1000])
colorbar

%% Calculate Pearson correlation coefficient between each component and the signal
MI = zeros(1, size(u,1));
for i = 1:size(u,1)
    MI(i) = corr(u(i,:)', y', 'type', 'Pearson');
end
disp(['Pearson correlation coefficient: ', num2str(MI)])

%% Correlation coefficient greater than threshold Thr
Thr = median(MI);
[mm, nn] = find(MI > Thr);

%% Denoised signal
rs = sum(u(nn,:));

% %% Calculate evaluation metrics
% [snr1, rmse1] = EvaMetrix(y, rs);
% disp(['EMD SNR: ', num2str(snr1)])
% disp(['EMD RMSE: ', num2str(rmse1)])
% CR_emds = corr(y', rs', 'type', 'Pearson');
% disp(['EMD Corr: ', num2str(CR_emds)])

% %% Plot curve
% figure
% plot(t, y, 'b-')
% hold on
% plot(t, ys, 'r-.')
% hold on
% plot(t, rs, 'g-s')
% xlabel('Time')
% ylabel('Amplitude')
% axis tight
% legend('Original Signal', 'Noisy Signal', 'EMD Denoised Signal')

end
