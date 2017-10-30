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
window_wavelength = 6;
overlap_ratio = 0.5;

window_size = round(window_wavelength * (1 / fc) * fs);

window = Windows(pre_full, post_full, window_size, overlap_ratio);
final = false;
t = 0;
delay = cell(1, round(length(pre_full) / window_size));
while ~final
%while t < 30;
    [pre, post, final] = window.Next();
    [val, lag] = xcorr(pre, post);
    [~, idx] = max(val);
    delay{t+1} = lag(idx);
    %{
    if t == 30
        plot(lag, val)
        break
    end
    %}
    t = t + 1;
end
delay

