clear;
data = load('FUS_RFData.mat');

pre_full = data.FUS_pre;
post_full = data.FUS_post;
fs = data.fs;
fc = data.fc;
window_wavelength = 10;
overlap_ratio = 0.75;

window_size = round(window_wavelength * (1 / fc) * fs);

window = Windows(pre_full, post_full, window_size, overlap_ratio);
final = false;
t = 0;

while ~final
%while t < 30;
    [pre, post, final] = window.Next();
end
