%
% EE6265 Fu-En Wang 106061531 HW2 11/14/2017
%

clear;
close all;
N_lst = [500 1000 2000 5000 10000];
M_lst = [20 10 5 2 1];
path = '../doc/src/2pi';
fig_flag = true;
bins = 20;

origin_amp = rand(1, 10000);
origin_phase = rand(1, 10000) * 2 * pi;
origin_data = origin_amp .* exp(1i * origin_phase);
origin_intensity = origin_amp.^2;

if fig_flag
    fig = figure();
    subplot('211')
    histogram(origin_amp, bins);
    tmp = sprintf('Amplitude Histogram(origin)');
    title(tmp)
    xlabel('value')
    ylabel('histogram')
    subplot('212')
    histogram(origin_intensity, bins);
    tmp = sprintf('Intensity Histogram(origin)');
    title(tmp)
    xlabel('value')
    ylabel('histogram')
    fig_name = sprintf('%s/hist_origin.pdf', path);
    saveFig(fig, fig_name);
end

all_mean_amp = zeros(1, length(M_lst));
all_std_amp = zeros(1, length(M_lst));
all_mean_intensity = zeros(1, length(M_lst));
all_std_intensity = zeros(1, length(M_lst));
ratio_amp = zeros(1, length(M_lst));
ratio_intensity = zeros(1, length(M_lst));

all_mean_amp_smooth_first = zeros(1, length(M_lst));
all_std_amp_smooth_first = zeros(1, length(M_lst));
all_mean_intensity_smooth_first = zeros(1, length(M_lst));
all_std_intensity_smooth_first = zeros(1, length(M_lst));
ratio_amp_smooth_first = zeros(1, length(M_lst));
ratio_intensity_smooth_first = zeros(1, length(M_lst));

all_mean_amp_smooth_later = zeros(1, length(M_lst));
all_std_amp_smooth_later = zeros(1, length(M_lst));
all_mean_intensity_smooth_later = zeros(1, length(M_lst));
all_std_intensity_smooth_later = zeros(1, length(M_lst));
ratio_amp_smooth_later = zeros(1, length(M_lst));
ratio_intensity_smooth_later = zeros(1, length(M_lst));

F = [0.5 0.5];

for i = 1:length(N_lst)
    N = N_lst(i);
    N_s = sprintf('%d', N_lst(i));
    
    M = M_lst(i);
    new_data = getNewArray(origin_data, M, N);
    new_data_smooth = filter(F, 1, new_data);
    new_amp = abs(new_data);
    new_intensity = new_amp.^2;
    
    new_amp_smooth_first = abs(new_data_smooth);
    new_intensity_smooth_first = new_amp_smooth_first.^2;
    
    new_amp_smooth_later = filter(F, 1,new_amp);
    new_intensity_smooth_later = filter(F, 1, new_intensity);
    
    all_mean_amp(i) = mean(new_amp);
    all_std_amp(i) = std(new_amp);
    all_mean_intensity(i) = mean(new_intensity);
    all_std_intensity(i) = std(new_intensity);
    
    all_mean_amp_smooth_first(i) = mean(new_amp_smooth_first);
    all_std_amp_smooth_first(i) = std(new_amp_smooth_first);
    all_mean_intensity_smooth_first(i) = mean(new_intensity_smooth_first);
    all_std_intensity_smooth_first(i) = std(new_intensity_smooth_first);
    
    all_mean_amp_smooth_later(i) = mean(new_amp_smooth_later);
    all_std_amp_smooth_later(i) = std(new_amp_smooth_later);
    all_mean_intensity_smooth_later(i) = mean(new_intensity_smooth_later);
    all_std_intensity_smooth_later(i) = std(new_intensity_smooth_later);
    
    if fig_flag
        fig = figure();
        subplot('211')
        histogram(new_amp, bins);
        tmp = sprintf('Amplitude Histogram(N = %d, M = %d)', N, M);
        title(tmp)
        xlabel('value')
        ylabel('histogram')
        subplot('212')
        histogram(new_intensity, bins);
        tmp = sprintf('Intensity Histogram(N = %d, M = %d)', N, M);
        title(tmp)
        xlabel('value')
        ylabel('histogram')
        fig_name = sprintf('%s/hist_%d_%d.pdf', path, N, M);
        saveFig(fig, fig_name);
    end
    
end
ratio_amp(:) = all_mean_amp ./ all_std_amp;
ratio_intensity(:) = all_mean_intensity ./ all_std_intensity;

ratio_amp_smooth_first(:) = all_mean_amp_smooth_first ./ all_std_amp_smooth_first;
ratio_intensity_smooth_first(:) = all_mean_intensity_smooth_first ./ all_std_intensity_smooth_first;

ratio_amp_smooth_later(:) = all_mean_amp_smooth_later ./ all_std_amp_smooth_later;
ratio_intensity_smooth_later(:) = all_mean_intensity_smooth_later ./ all_std_intensity_smooth_later;

% plot ratio
if fig_flag
    fig = figure();
    subplot('211');
    fig_title = 'Ratio (Amplitude)';
    plot(M_lst, ratio_amp);
    title(fig_title)
    xlabel('M')
    ylabel('ratio')
    subplot('212');
    fig_title = 'Ratio (Intensity)';
    plot(M_lst, ratio_intensity);
    title(fig_title)
    xlabel('M')
    ylabel('ratio')
    fig_path = [path '/ratio'  '.pdf'];
    saveFig(fig, fig_path);
end
% plot smooth_later for problem (d) 2 x 1 subplot
if fig_flag
    fig = figure(); 
    subplot(2, 1, 1);
    fig_title = 'Ratio (Smoothed Amplitude)';
    plot(M_lst, ratio_amp_smooth_later);
    title(fig_title)
    xlabel('M')
    ylabel('ratio')
    
    subplot(2, 1, 2);
    fig_title = 'Ratio (Smoothed Intensity)';
    plot(M_lst, ratio_amp_smooth_later);
    title(fig_title)
    xlabel('M')
    ylabel('ratio')
    name = sprintf('%s/ratio_smooth_later.pdf', path);
    saveFig(fig, name);
end
% plot smooth_first for problem (e) 5 x 2 subplot
if fig_flag
    fig = figure();
    
    subplot(2, 1, 1);
    fig_title = 'Ratio (Amplitude of smoothed data)';
    plot(M_lst, ratio_amp_smooth_first);
    title(fig_title)
    xlabel('M')
    ylabel('ratio')
    
    subplot(2, 1, 2);
    fig_title = 'Ratio (Intensity of smoothed data)';
    plot(M_lst, ratio_amp_smooth_first);
    title(fig_title)
    xlabel('M')
    ylabel('ratio')
    
    name = sprintf('%s/ratio_smooth_first.pdf', path);
    saveFig(fig, name);
end

