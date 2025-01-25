function rs=VMDMain(time,fs,y,ns)
%% VMD Parameter Settings
alpha=2000;
K=8;
tau = 0;           
DC = 0;             
init = 1;          
tol = 1e-7;

%% VMD
[u, u_hat, omega] = VMD(ns, alpha, tau, K, DC, init, tol);
[m,~]=size(u);

% figure
% for i=1:m-1
% subplot(m,1,i)
% plot(time,u(i,:),'b-','linewidth',1)
% ylabel(['IMF',num2str(i)]);
% axis tight
% end
% subplot(m,1,m)
% plot(time,u(m,:),'b-','linewidth',1)
% ylabel('Res');
% axis tight
% 
% figure
% for i=1:m-1
% subplot(m,1,i)

% %% FFT
% [cc,y_f]=hua_fft(u(i,:),fs,1);
% plot(y_f,cc,'b','LineWIdth',1.5);
% ylabel(['FFT IMF',num2str(i)]);
% axis tight
% end
% subplot(m,1,m)
% [cc,y_f]=hua_fft(u(m,:),fs,1);
% plot(y_f,cc,'b','LineWIdth',1.5);
% axis tight
% ylabel('FFT Res');

%% Calculate Envelope Entropy
[m,n]=size(u);
for i=1:m
    En(i)=BaoluoEntropy(u(i,:));
end
disp(['Envelope entropy of each IMF component: ' num2str(En)])

%% Calculate the Pearson coefficient between each component and the signal
MI=zeros(1,size(u,1));
for i = 1:size(u,1)
    MI(i) = corr(u(i,:)',y','type','Pearson');
end
disp(['Pearson correlation coefficients between VMD components and the signal: ', num2str(MI)])

%% Correlation coefficient greater than threshold Thr
Thr=median(MI);
[~,nn]=find(MI>Thr);
%% Denoised Signal
rs=sum(u(nn,:));
% %% Calculate Evaluation Metrics
% [snr1,rmse1]=EvaMetrix(y,rs);
% CR=corr(y',rs','type','Pearson');
% disp(['VMD SNR: ' num2str(snr1)])
% disp(['VMD RMSE: ' num2str(rmse1)])
% disp(['VMD Corr: ' num2str(CR)])
