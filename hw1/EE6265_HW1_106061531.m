%
% EE6265 ??? 106061531 HW1 10/30/2017
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
window_wavelength = 10;
overlap_ratio = 0;

window_size = round(window_wavelength * (1 / fc) * fs);

window = Windows(pre_full, post_full, window_size, overlap_ratio);
final = false;
t = 0;
delay = zeros(1, round(length(pre_full) / window_size));
center_idx = zeros(1, round(length(pre_full) / window_size));
while ~final
%while t < 30;
    [pre, post, center, final] = window.Next();
    [val, lag] = xcorr(pre, post);
    [~, idx] = max(val);
    delay(t+1) = lag(idx);
    center_idx(t+1) = center;
    %{
    if t == 30
        plot(lag, val)
        break
    end
    %}
    t = t + 1;
end

depth = c0 * center_idx * (1 / fs);
delay_sec = delay * (1 / fs);
tmp1 = [0 depth(1:end-1)];
tmp2 = [0 delay_sec(1:end-1)];

delta_depth = depth - tmp1;
delta_delay = delay_sec - tmp2;

figure()
plot(depth * 1e3, delta_delay ./ delta_depth);
title('Strain')
xlabel('Depth(mm)')
ylabel('Strain')
% figure()
% title('$a = \frac{4}{3}$', 'interpreter', 'latex')

