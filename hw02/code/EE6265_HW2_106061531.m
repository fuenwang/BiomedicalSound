%
% EE6265 Fu-En Wang 106061531 HW2 11/14/2017
%
clear;
close all;
N_lst = [500 1000 2000 5000 10000];
M_lst = [1 2 5 10 20];
path = '../doc/src/2pi';
fig_flag = true;

origin_amp = rand(1, 10000);
origin_phase = rand(1, 10000) * 2 * pi;
origin_data = origin_amp .* exp(1i * origin_phase);
origin_intensity = origin_amp.^2;

if fig_flag
    fig = figure();
    subplot('211')
    histogram(origin_amp);
    tmp = sprintf('Amplitude Histogram(origin)');
    title(tmp)
    xlabel('value')
    ylabel('histogram')
    subplot('212')
    histogram(origin_intensity);
    tmp = sprintf('Intensity Histogram(origin)');
    title(tmp)
    xlabel('value')
    ylabel('histogram')
    fig_name = sprintf('%s/hist_origin.pdf', path);
    saveFig(fig, fig_name);
end

all_mean_amp = zeros(length(N_lst), length(M_lst));
all_std_amp = zeros(length(N_lst), length(M_lst));
all_mean_intensity = zeros(length(N_lst), length(M_lst));
all_std_intensity = zeros(length(N_lst), length(M_lst));
ratio_amp = zeros(length(N_lst), length(M_lst));
ratio_intensity = zeros(length(N_lst), length(M_lst));

all_mean_amp_smooth_first = zeros(length(N_lst), length(M_lst));
all_std_amp_smooth_first = zeros(length(N_lst), length(M_lst));
all_mean_intensity_smooth_first = zeros(length(N_lst), length(M_lst));
all_std_intensity_smooth_first = zeros(length(N_lst), length(M_lst));
ratio_amp_smooth_first = zeros(length(N_lst), length(M_lst));
ratio_intensity_smooth_first = zeros(length(N_lst), length(M_lst));

all_mean_amp_smooth_later = zeros(length(N_lst), length(M_lst));
all_std_amp_smooth_later = zeros(length(N_lst), length(M_lst));
all_mean_intensity_smooth_later = zeros(length(N_lst), length(M_lst));
all_std_intensity_smooth_later = zeros(length(N_lst), length(M_lst));
ratio_amp_smooth_later = zeros(length(N_lst), length(M_lst));
ratio_intensity_smooth_later = zeros(length(N_lst), length(M_lst));

F = [0.5 0.5];

for i = 1:length(N_lst)
    N = N_lst(i);
    N_s = sprintf('%d', N_lst(i));
    for j = 1:length(M_lst)
        M = M_lst(j);
        new_data = getNewArray(origin_data, M, N);
        new_data_smooth = filter(F, 1, new_data);
        new_amp = abs(new_data);
        new_intensity = new_amp.^2;
        
        new_amp_smooth_first = abs(new_data_smooth);
        new_intensity_smooth_first = new_amp_smooth_first.^2;
        
        new_amp_smooth_later = filter(F, 1,new_amp);
        new_intensity_smooth_later = filter(F, 1, new_intensity);
        
        all_mean_amp(i, j) = mean(new_amp);
        all_std_amp(i, j) = std(new_amp);
        all_mean_intensity(i, j) = mean(new_intensity);
        all_std_intensity(i, j) = std(new_intensity);
        
        all_mean_amp_smooth_first(i, j) = mean(new_amp_smooth_first);
        all_std_amp_smooth_first(i, j) = std(new_amp_smooth_first);
        all_mean_intensity_smooth_first(i, j) = mean(new_intensity_smooth_first);
        all_std_intensity_smooth_first(i, j) = std(new_intensity_smooth_first);
        
        all_mean_amp_smooth_later(i, j) = mean(new_amp_smooth_later);
        all_std_amp_smooth_later(i, j) = std(new_amp_smooth_later);
        all_mean_intensity_smooth_later(i, j) = mean(new_intensity_smooth_later);
        all_std_intensity_smooth_later(i, j) = std(new_intensity_smooth_later);
        
        if (M == 1 || M == 10) && fig_flag
            fig = figure();
            subplot('211')
            histogram(new_amp);
            tmp = sprintf('Amplitude Histogram(N = %d, M = %d)', N, M);
            title(tmp)
            xlabel('value')
            ylabel('histogram')
            subplot('212')
            histogram(new_intensity);
            tmp = sprintf('Intensity Histogram(N = %d, M = %d)', N, M);
            title(tmp)
            xlabel('value')
            ylabel('histogram')
            fig_name = sprintf('%s/hist_%d_%d.pdf', path, N, M);
            saveFig(fig, fig_name);
        end
    end
end
ratio_amp(:, :) = all_mean_amp ./ all_std_amp;
ratio_intensity(:, :) = all_mean_intensity ./ all_std_intensity;

ratio_amp_smooth_first(:, :) = all_mean_amp_smooth_first ./ all_std_amp_smooth_first;
ratio_intensity_smooth_first(:, :) = all_mean_intensity_smooth_first ./ all_std_intensity_smooth_first;

ratio_amp_smooth_later(:, :) = all_mean_amp_smooth_later ./ all_std_amp_smooth_later;
ratio_intensity_smooth_later(:, :) = all_mean_intensity_smooth_later ./ all_std_intensity_smooth_later;

% plot ratio
for i = 1:length(N_lst)
    if fig_flag
        N_s = sprintf('%d', N_lst(i));
        fig = figure();
        subplot('211');
        fig_title = sprintf('N = %d (Amplitude)', N_lst(i));
        plot(M_lst, ratio_amp(i, :));
        title(fig_title)
        xlabel('M')
        ylabel('ratio')
        subplot('212');
        fig_title = sprintf('N = %d (Intensity)', N_lst(i));
        plot(M_lst, ratio_intensity(i, :));
        title(fig_title)
        xlabel('M')
        ylabel('ratio')
        fig_path = [path '/ratio_' N_s '.pdf'];
        saveFig(fig, fig_path);
    end
end

% plot smooth_later for problem (d) 5 x 2 subplot
if fig_flag
    fig = figure();
    for i = 1:length(N_lst)
        left_id = i * 2 - 1;
        right_id = left_id + 1;

        subplot(5, 2, left_id);
        fig_title = sprintf('N = %d (Smoothed Amplitude)', N_lst(i));
        plot(M_lst, ratio_amp_smooth_later(i, :));
        title(fig_title)
        xlabel('M')
        ylabel('ratio')
        
        subplot(5, 2, right_id);
        fig_title = sprintf('N = %d (Smoothed Intensity)', N_lst(i));
        plot(M_lst, ratio_amp_smooth_later(i, :));
        title(fig_title)
        xlabel('M')
        ylabel('ratio')
    end
    set(fig, 'Position', [0 0 800 1080]);
    name = sprintf('%s/ratio_smooth_later.pdf', path);
    saveFig(fig, name);
end

% plot smooth_first for problem (e) 5 x 2 subplot
if fig_flag
    fig = figure();
    for i = 1:length(N_lst)
        left_id = i * 2 - 1;
        right_id = left_id + 1;

        subplot(5, 2, left_id);
        fig_title = sprintf('N = %d (Amplitude of smoothed data)', N_lst(i));
        plot(M_lst, ratio_amp_smooth_first(i, :));
        title(fig_title)
        xlabel('M')
        ylabel('ratio')
        
        subplot(5, 2, right_id);
        fig_title = sprintf('N = %d (Intensity of smoothed data)', N_lst(i));
        plot(M_lst, ratio_amp_smooth_first(i, :));
        title(fig_title)
        xlabel('M')
        ylabel('ratio')
    end
    set(fig, 'Position', [0 0 800 1080]);
    name = sprintf('%s/ratio_smooth_first.pdf', path);
    saveFig(fig, name);
end

