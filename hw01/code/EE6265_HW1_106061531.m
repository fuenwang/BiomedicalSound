%
% EE6265 Fu-En Wang 106061531 HW1 10/30/2017
%
% This file will run the flow to estimate time-shift of signal.
% There is another file "Windows.m", it will create an object
% to mantain windows.
%

%clear;

data = load('FUS_RFData.mat');

pre_full = data.FUS_pre;
post_full = data.FUS_post;
fs = data.fs * 1e6;
fc = data.fc * 1e6;
c0 = data.c0 * 1e-3 / 1e-6;

window_wavelength = 2;
overlap_ratio = 0.75;

resample_factor = 10;

pre_full = interp(pre_full, resample_factor);
post_full = interp(post_full, resample_factor);
fs = fs * resample_factor;
window_size = round(window_wavelength * (1 / fc) * fs);

window = Windows(pre_full, post_full, window_size, overlap_ratio);
final = false;
t = 0;
delay = zeros(1, round(length(pre_full) / window_size));
center_idx = zeros(1, round(length(pre_full) / window_size));
while ~final
    [pre, post, center, final] = window.Next();
%     if t == 1
%         pre(end - round(window_size * overlap_ratio): end)
%     elseif t == 2
%         pre(1:round(window_size * overlap_ratio))
%    end
    [val, lag] = xcorr(pre, post);
    [~, idx] = max(val);
    delay(t+1) = lag(idx);
    center_idx(t+1) = center;
    t = t + 1;
end
new_fs = length(delay) / length(pre_full) * fs;
cutoff =  4.398 * 1e4;
avg_filter_size = round(sqrt((0.442947 * new_fs / cutoff)^2 + 1));
%avg_filter_size = 30;
delay_sec = delay * (1 / fs);
delay_sec = filter(ones(1, avg_filter_size) / avg_filter_size, 1, delay_sec);
depth = center_idx * (1 / fs) / 2 * c0;
strain = diff(delay_sec * c0 / 2) ./ diff(depth);

% x = fft(delay_sec - mean(delay_sec));
% x = fftshift(x);
% side = (length(x) - 1) / 2;
% f_axis = (-side:side) / 2 / side * new_fs;
% fig = figure();
% plot(f_axis, abs(x), 'LineWidth', 1);
% title('FFT of Echo-Time Shift(smoothed)')
% xlabel('Hz')
% ylabel('Magnitude')
% fig.PaperPositionMode = 'auto';
% fig_pos = fig.PaperPosition;
% fig.PaperSize = [fig_pos(3) fig_pos(4)];
% print(fig, '../doc/src/fft.pdf', '-dpdf')
% error

figure()
plot(depth*1e3, (delay_sec / 1) * 1e6)
title('Echo-Time Shift')
xlabel('Depth(mm)')
ylabel('Delay(us)')
%error
figure()
plot(depth(1:end-1)*1e3, strain * 100)
title('Thermal Strain(%)')
xlabel('Depth(mm)', 'FontSize', 15)
ylabel('\boldmath{$\frac{\partial \Delta t(z)}{\partial z}$ (\%)}', 'interpreter', 'latex', 'FontSize', 20)








