clear;
data = load('FUS_RFData.mat');

pre_full = data.FUS_pre;
post_full = data.FUS_post;
fs = data.fs;
window_size = 10;


window = Windows(pre_full, post_full, window_size);
final = false;
while ~final
    [pre, post, final] = window.Next();
end