%
% EE6265 Fu-En Wang 106061531 HW1 10/30/2017
%
% This file will run the flow to estimate time-shift of signal.
% There is another file "Windows.m", it will create an object
% to mantain windows.
%

clear;
data = load('FUS_RFData.mat');

pre_full = data.FUS_pre;
post_full = data.FUS_post;
fs = data.fs * 1e6;
fc = data.fc * 1e6;
c0 = data.c0 * 1e-3 / 1e-6;

window_wavelength = 2;
overlap_ratio = 0.5;
resample_factor = 3;
avg_filter_size = 5;

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
    [val, lag] = xcorr(pre, post);
    [~, idx] = max(val);
    delay(t+1) = lag(idx);
    center_idx(t+1) = center;
    t = t + 1;
end
delay_sec = delay * (1 / fs);
delay_sec = filter(ones(1, avg_filter_size) / avg_filter_size, 1, delay_sec);
depth = center_idx * (1 / fs) / 2 * c0;

figure()
plot(depth*1e3, (delay_sec / 2) * 1e6)
title('Echo-Time Shift')
xlabel('Depth(mm)')
ylabel('Delay(us)')








