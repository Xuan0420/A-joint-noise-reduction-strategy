t = 0:0.1:100;
y = 1/10*t - 5;

% Threshold selection
thr = thselect(y, 'minimaxi');
% g = std(y)
% thr = sqrt(21*g*100);

% Calculate soft and hard threshold functions
sorh_s = 's'; % Select soft threshold
sorh_h = 'h'; % Select hard threshold
soft_thr = wthresh(y, sorh_s, thr);
hard_thr = wthresh(y, sorh_h, thr);
imd_thr_1 = improved_threshold2(y, thr, 1, 1);
imd_thr_2 = improved_threshold2(y, thr, 10, 0.1);
imd_thr_3 = improved_threshold2(y, thr, 1, 100);

% Plot traditional soft and hard threshold functions
figure;

subplot(1, 2, 1);
plot(y, hard_thr, 'LineWidth', 1);
title('Hard');
text(0.05, 0.9, '(a)', 'Units', 'Normalized', 'FontSize', 12);
xlabel('w');
ylabel('w1');
grid on; % Add grid lines

subplot(1, 2, 2);
plot(y, soft_thr, 'LineWidth', 1);
title('Soft');
text(0.05, 0.9, '(b)', 'Units', 'Normalized', 'FontSize', 12);
xlabel('w');
ylabel('w1');
grid on; % Add grid lines


% % Plot the original signal
figure;
plot(y, y, 'w', 'LineWidth', 1.5);
hold on; % Enable hold state for plotting

% Plot signal after hard thresholding
plot(y, hard_thr, 'color', [014/255 062/255 135/255], 'LineWidth', 1.5);
% Plot signal after soft thresholding
plot(y, soft_thr, 'color', [216/255 178/255 058/255], 'LineWidth', 1.5);
% Plot signal after improved thresholding
plot(y, imd_thr_3, 'color', [052/255 108/255 172/255], 'LineWidth', 1.5);
% Plot signal after improved thresholding
plot(y, imd_thr_2, 'color', [222/255 234/255 234/255], 'LineWidth', 1.5);
% Plot signal after improved thresholding
plot(y, imd_thr_1, 'color', [247/255 228/255 116/255], 'LineWidth', 1.5);

hold off; % Disable hold state for plotting

legend('', 'Hard', 'Soft', 'Improved(α=1,β=100)', 'Improved(α=10,β=0.1)', 'Improved(α=1,β=1)');
% title('Improved threshold functions');
xlabel('w');
ylabel('w1');
grid on;
