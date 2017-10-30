clear;
data = load('FUS_RFData.mat');

pre_full = data.FUS_pre;
post_full = data.FUS_post;
fs = data.fs;
fc = data.fc;
window_wavelength = 6;
overlap_ratio = 0;

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
