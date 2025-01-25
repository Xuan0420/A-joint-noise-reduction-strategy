%%%%%% SMA-VMD-IWTD
clc
clear 
close all
warning off
%% Parameter settings
fs = 4096;  % Sampling frequency
t_total = 0.2;  % Total time
SNR = 5;  % Signal-to-noise ratio of the measurement noise
t = 0:1/fs:t_total;  % Time range
%% System parameters
R = [2, 2, 2, 2];  % Amplitude
zeta = [0.3661, 0.0438, 0.0368, 0.0062];  % Damping coefficient
w_m = 2*pi*[49.23, 172.76, 260.03, 625.93];  % Angular frequency
omega_d = w_m.*sqrt(1-zeta.^2);  % Damped oscillation frequency
%% Excitation signal
f = exp(-t);  
%% Response signal
y = zeros(size(t));
for r = 1:length(R)
    y = y + R(r)*exp(-zeta(r)*w_m(r)*t).*sin(omega_d(r)*t);
end
%% Add Gaussian white noise
noise_power = var(y)/10^(SNR/10);
randn('state', 0);  % Set random seed for consistent results
%% Noise signal
w = sqrt(noise_power)*randn(size(t));
ys = y + w;
%% Calculate evaluation metrics
[snr0,psnr0,rmse0]=EvaMetrix(y,ys);
disp(['Signal-to-noise ratio (SNR) of the original noisy signal: ' num2str(snr0)])
disp(['Peak signal-to-noise ratio (PSNR) of the original noisy signal: ' num2str(psnr0)])
disp(['Root mean square error (RMSE) of the original noisy signal: ' num2str(rmse0)])

% %% Plot the curve
% figure
% plot(t,y,'b-','linewidth',1)
% hold on
% plot(t,ys,'r-','linewidth',1)
% xlabel('Time')
% ylabel('Amplitude')
% axis tight
% legend('Original signal','The signal with noise')

%% SMA
%% Parameter settings
pop=20; % Population size
Maxgen=20; % Maximum number of iterations
%% Parameter range (alpha and K)
ub=[2500 10];
lb=[100 3];
%% Optimization variable dimensions
dim =2;   % 2 unknowns
%% Fitness function
fobj=@(X)objfun(X,ys);
%% Slime Mould Algorithm
[SMA_worst_score, SMA_worst_pos, SMA_Curve] = SMA(pop, Maxgen, lb, ub, dim, fobj);
%% Convergence curve
figure
plot(SMA_Curve,'-*b','linewidth',1.5)
title('SMA convergence curve');
xlabel('Number of iterations')
ylabel(['Fitness'])
axis([1 length(SMA_Curve) 5.850 5.870]);
%axis tight
%% Optimal parameters
xxx=SMA_worst_pos;
%% Assign optimal parameters to alpha and K
alpha=round(xxx(1));
K=round(xxx(2));
disp(alpha)
disp(K)
tau = 0;           
DC = 0;             
init = 1;          
tol = 1e-7;
%% SMA-VMD
[u, u_hat, omega] = VMD(ys, alpha, tau, K, DC, init, tol);
[m,~]=size(u);
%% Plot SMA-VMD decomposition
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
%% FFT transform
[cc,y_f]=hua_fft(u(i,:),fs,1);
plot(y_f,cc,'b','LineWidth',1.5);
ylabel(['FFT IMF',num2str(i)]);
axis tight
end
subplot(m,1,m)
[cc,y_f]=hua_fft(u(m,:),fs,1);
plot(y_f,cc,'b','LineWidth',1.5);
axis tight
ylabel('FFT Res');

figure;
hold on;
for i = 1:m
    [cc,y_f]=hua_fft(u(i,:),fs,1);
    plot(y_f, cc, 'LineWidth', 1);
axis tight;
end
hold off;
xlabel('Frequency');
ylabel('Amplitude');
legend('IMF 1', 'IMF 2', 'IMF 3','IMF 4','IMF 5','IMF 6','IMF 7','IMF 8','IMF 9','Res');

figure;
hold on;
for i = 1:m
    [cc,y_f]=hua_fft(u(i,:),fs,1);
    plot3(i*ones(size(y_f)), y_f, cc);
end
hold off;
xlabel('IMF');
ylabel('Frequency');
zlabel('Amplitude');
legend('IMF 1', 'IMF 2', 'IMF 3','IMF 4','IMF 5','IMF 6','IMF 7','IMF 8','IMF 9','Res');

%% Calculate envelope entropy
[m,n]=size(u);
for i=1:m
    En(i)=BaoluoEntropy(u(i,:));
end
disp(['Envelope entropy of each IMF component: ' num2str(En)])

%% Calculate Pearson correlation coefficient between each component and the signal
MI=zeros(1,size(u,1));
for i = 1:size(u,1)
    MI(i) = corr(u(i,:)',y','type','Pearson');
end
disp(['Pearson correlation coefficients between SMA-VMD components and the signal: ',num2str(MI)])
%% Correlation coefficient greater than threshold Thr 
Thr=median(MI);
[mm,nn]=find(MI>Thr);

%% Calculate power spectral entropy of each component
pse = psd_entropy(u,fs);
disp(pse);

% figure;
% plot(pse, 'or-', 'LineWidth', 1.5);
% axis([1 length(psd) 0 1]);
% xlabel('IMF');
% ylabel('PSE');

figure;
hold on;
plot(MI, 'ok-', 'MarkerSize', 6, 'LineWidth', 1.5);
plot(pse, 'or-', 'MarkerSize', 6, 'LineWidth', 1.5);
hold off;
xlabel('IMF');
legend('PC','PSE');
axis tight;
xlim([1 length(pse)]);
ylim([0 1]);
box on;

%% HHT spectrogram
[A,f,tt] = hhspectrum(u); % Apply Hilbert transform to each IMF component and then perform FFT; f is normalized frequency
%% A is instantaneous amplitude, f is instantaneous frequency
%% Combine HHT spectrogram
[E,ttt,Cenf]=toimage(A,f);
%% Plot HHT spectrogram
figure
imagesc(t,[0,0.5*fs],E);
set(gca,'YDir','normal')
%text(0.05, 0.9, '(b)', 'Units', 'Normalized', 'FontSize', 14,'Color', 'white');
xlabel('t/s')
ylabel('f/Hz')
title('Hilbert spectrum (VMD)')
ylim([0,1000])
colorbar

% %% Black and white HHT spectrogram
% figure
% imagesc(t,[0,0.5*fs],E);
% set(gca,'YDir','normal')
% xlabel('time')
% ylabel('frequency/Hz')
% title('Hilbert spectrum')
% ylim([0,1000])
% colormap(flipud(hot)); 
% colorbar


%% SMA-VMD denoised signal
rs=sum(u(nn,:));
%% Calculate evaluation metrics
[snr1,psnr1,rmse1]=EvaMetrix(y,rs);
disp(['SMA-VMD SNR: ' num2str(snr1)])
disp(['SMA-VMD PSNR: ' num2str(psnr1)])
disp(['SMA-VMD RMSE: ' num2str(rmse1)])

%% SMA-VMD-WTD (Soft Threshold)
wavelet = 'db4'; % Wavelet basis function
level = 8; % Decomposition level
threshold = 0.08; % 10db, thr=0.08
rs_ctjb = softThresholdDenoise(rs, wavelet, level, threshold);

%% Calculate evaluation metrics
[snr2,psnr2,rmse2]=EvaMetrix(y,rs_ctjb);
disp(['SMA-VMD-WTD SNR: ' num2str(snr2)])
disp(['SMA-VMD-WTD PSNR: ' num2str(psnr2)])
disp(['SMA-VMD-WTD RMSE: ' num2str(rmse2)])

%% SMA-VMD-IWTD (Improved Threshold)
alpha = 10; 
beta = 0.1; 
rs_gjxb = improvedWaveletDenoise(rs, wavelet, level, threshold, alpha, beta);
%% Calculate evaluation metrics
[snr3,psnr3,rmse3]=EvaMetrix(y,rs_gjxb);
disp(['SMA-VMD-IWTD SNR: ' num2str(snr3)])
disp(['SMA-VMD-IWTD PSNR: ' num2str(psnr3)])
disp(['SMA-VMD-IWTD RMSE: ' num2str(rmse3)])


%% Ordinary VMD
rrs=VMDMain(t,fs,y,ys);
%% Calculate evaluation metrics
[snr4,psnr4,rmse4]=EvaMetrix(y,rrs);
disp(['VMD SNR: ' num2str(snr4)])
disp(['VMD PSNR: ' num2str(psnr4)])
disp(['VMD RMSE: ' num2str(rmse4)])

%% Ordinary WTD
wtds = softThresholdDenoise(ys, wavelet, level, threshold);
%% Calculate evaluation metrics
[snr5,psnr5,rmse5]=EvaMetrix(y,wtds);
disp(['WTD SNR: ' num2str(snr5)])
disp(['WTD PSNR: ' num2str(psnr5)])
disp(['WTD RMSE: ' num2str(rmse5)])


%% Ordinary EMD
emds=EMDMain(t,fs,y,ys);
%% Calculate evaluation metrics
[snr6,psnr6,rmse6] = EvaMetrix(y,emds);
disp(['EMD SNR: ' num2str(snr6)])
disp(['EMD PSNR: ' num2str(psnr6)])
disp(['EMD RMSE: ' num2str(rmse6)])

%% Plot method effect curve
%% Plotting method effect curves
figure
% plot(t,ys,'b','LineWIdth',1) %ys is the noisy signal
% hold on
plot(t,emds,'color', [030/255 070/255 110/255],'LineWIdth',1) %emds is the signal after EMD denoising
hold on
plot(t,wtds,'color', [114/255 188/255 213/255],'LineWIdth',1) %wtds is the signal after WTD denoising
hold on
plot(t,rrs,'color', [170/255 220/255 224/255],'LineWIdth',1) %rrs is the signal after VMD denoising
hold on
plot(t,rs,'color', [225/255 230/255 183/255],'LineWIdth',1) %rs is the signal after SMA-VMD denoising
hold on
plot(t,rs_ctjb,'y','LineWIdth',1) %rs_ctjb is the signal after SMA-VMD-WTD denoising
hold on
plot(t,rs_gjxb,'color', [247/255 170/255 088/255],'LineWIdth',1) %rs_gjxb is the signal after SMA-VMD-IWTD denoising
hold on
plot(t,y,'r--','LineWIdth',1) %y is the original signal
hold on

xlabel('t')
ylabel('A')
axis tight
legend('EMD denoising signal','WTD denoising signal','VMD denoising signal','SMA-VMD denoising signal','SMA-VMD-WTD denoising signal','SMA-VMD-IWTD denoising signal','Original signal');
title('Denoising effects of methods')


% % Plot the first subplot, original signal vs noisy signal
% subplot(3, 3, 1);
% plot(t, ys, 'b', 'LineWidth', 1);  % ys is the noisy signal
% hold on;
% plot(t, y, 'r--', 'LineWidth', 1);  % y is the original signal
% xlabel('Time');
% ylabel('Amplitude');
% title('Original Signal vs Noisy Signal');

% Plot the second subplot, signal after EMD denoising
subplot(2, 3, 1);
plot(t, emds, 'Color', [30/255, 70/255, 110/255], 'LineWidth', 1);  % emds is the signal after EMD denoising
hold on;
plot(t, y, 'm--', 'LineWidth', 1);  % y is the original signal
xlabel('t');
ylabel('A');
legend('EMD','Ideal signal')
title('EMD');

% Plot the third subplot, signal after WTD denoising
subplot(2, 3, 2);
plot(t, wtds, 'Color', [30/255, 70/255, 110/255], 'LineWidth', 1);  % wtds is the signal after WTD denoising
hold on;
plot(t, y, 'm--', 'LineWidth', 1);  % y is the original signal
xlabel('t');
ylabel('A');
legend('WTD','Original')
title('WTD');

% Plot the fourth subplot, signal after VMD denoising
subplot(2, 3, 3);
plot(t, rrs, 'Color', [30/255, 70/255, 110/255], 'LineWidth', 1);  % rrs is the signal after VMD denoising
hold on;
plot(t, y, 'm--', 'LineWidth', 1);  % y is the original signal
xlabel('t');
ylabel('A');
legend('VMD','Ideal signal')
title('VMD');

% Plot the fifth subplot, signal after SMA-VMD denoising
subplot(2, 3, 4);
plot(t, rs, 'Color', [30/255, 70/255, 110/255], 'LineWidth', 1);  % rs is the signal after SMA-VMD denoising
hold on;
plot(t, y, 'm--', 'LineWidth', 1);  % y is the original signal
xlabel('t');
ylabel('A');
legend('SMA-VMD','Original')
title('SMA-VMD');

% Plot the sixth subplot, signal after SMA-VMD-WTD denoising
subplot(2, 3, 5);
plot(t, rs_ctjb,'Color',[30/255, 70/255, 110/255], 'LineWidth', 1);  % rs_ctjb is the signal after SMA-VMD-WTD denoising
hold on;
plot(t, y, 'm--', 'LineWidth', 1);  % y is the original signal
xlabel('t');
ylabel('A');
legend('SMA-VMD-WTD','Original')
title('SMA-VMD-WTD');

% Plot the seventh subplot, signal after SMA-VMD-IWTD denoising
subplot(2, 3, 6);
plot(t, rs_gjxb, 'Color', [30/255, 70/255, 110/255], 'LineWidth', 1);  % rs_gjxb is the signal after SMA-VMD-IWTD denoising
hold on;
plot(t, y, 'm--', 'LineWidth', 1);  % y is the original signal
legend('SMA-VMD-IWTD','Original')
xlabel('t');
ylabel('A');
title('SMA-VMD-IWTD');

% %% Plot method effect curve
% figure
% plot(t,y,'k','LineWIdth',1) % y is the original signal
% hold on
% % plot(t,ys,'b','LineWIdth',1) % ys is the noisy signal
% % hold on
% plot(t,emds,'LineWIdth',1) % emds is the signal after EMD denoising
% hold on
% plot(t,wtds,'LineWIdth',1) % wtds is the signal after WTD denoising
% hold on
% plot(t,rrs,'LineWIdth',1) % rrs is the signal after VMD denoising
% hold on
% plot(t,rs,'LineWIdth',1) % rs is the signal after SMA-VMD denoising
% hold on
% plot(t,rs_ctjb,'LineWIdth',1) % rs_ctjb is the signal after SMA-VMD-WTD denoising
% hold on
% plot(t,rs_gjxb, 'LineWIdth',1) % rs_gjxb is the signal after SMA-VMD-IWTD denoising
% hold on
% xlabel('Time')
% ylabel('Amplitude')
% axis tight
% legend('Original signal','EMD denoising signal','WTD denoising signal','VMD denoising signal','SMA-VMD denoising signal','SMA-VMD-WTD denoising signal','SMA-VMD-IWTD denoising signal');
% title('Denoising effects of methods')