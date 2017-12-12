%
%
%   EE6265 Fu-En Wang 106061531 HW3 12/10/2017
%
%
close all
clear;

fig_path = './';
GAMMA = 0.12;
H0 = 0.09; % J/cm^2
ua = 180; % cm^-1
speed = 1.5; % mm/us
thick = 300; % um
options = optimset('Display','off');

p1a = true;
p1b = true;
p1d = true;
p1e = true;
p1f = true;
p2a = true;
p2b = false;


if p2a && p2b
    error('Only one of p2a and p2b can be true')
end
%{
p1a = true;
p1b = true;
p1d = true;
p1e = true;
p1f = true;
%}
%
% Problem 1-a
%
time = -400:(thick/2); % ns
pos_pressure = ua * GAMMA * H0 * exp(-ua * 10^2 * abs(time) * speed * 10^-6);
neg_pressure = -flip(pos_pressure);
pressure = [pos_pressure neg_pressure];
[max_val, max_idx] = max(pressure);
[min_val, min_idx] = min(pressure);
time_axis = 0:length(pressure)-1;
A = [time_axis(max_idx) 1; time_axis(min_idx) 1];
b = [max_val min_val]';
x = A \ b;
pressure(min_idx:-1:max_idx) = time_axis(min_idx:-1:max_idx) * x(1) + x(2);
if p1a
    fig = figure();
    plot(time_axis, pressure, 'linewidth', 2)
    title('Acoustic wave')
    xlabel('time(ns)')
    ylabel('Pressure')
    saveFig(fig, [fig_path 'p1a.pdf'])
end

%
% Problem 1-b and Problem 1-c
%
gauss_std = (max_val) * 5 / 100;
p1a_time = time_axis;
p1a_pressure = pressure + gauss_std * randn(1, length(pressure));
time = 0:300;
pressure = ua * GAMMA * H0 * exp(-ua * 10^2 * abs(time) * speed * 10^-6);
pressure = pressure + gauss_std * randn(1, length(pressure));


x0 = 185;
func = @(x)sum(abs(pressure - x * GAMMA * H0 * exp(-x * abs(time) * speed * 10^-4)));
x = lsqnonlin(func,x0, [], [], options);
paper_pressure = max_val*exp(-185 * abs(time) * speed * 10^-4);
predict_pressure = max_val*exp(-x * abs(time) * speed * 10^-4);
if p1b
    fig = figure();
    plot(p1a_time, p1a_pressure, 'linewidth', 1)
    title('Acoustic wave')
    xlabel('time(ns)')
    ylabel('Pressure')
    saveFig(fig, [fig_path 'p1b-1.pdf'])
    
    fig = figure();
    plot(time, pressure, 'linewidth', 1)
    hold on
    plot(time, paper_pressure, 'r', 'linewidth', 2)
    legend('Signal', sprintf('%.2fexp(-%.2fz)', max_val, x0))
    title('Acoustic wave')
    xlabel('time(ns)')
    ylabel('Pressure')
    ylim([-0.3 2.1])
    saveFig(fig, [fig_path 'p1b-2.pdf'])
    
    fig = figure();
    plot(time, pressure, 'linewidth', 1)
    hold on
    plot(time, predict_pressure, 'r', 'linewidth', 2)
    legend('Signal', sprintf('%.2fexp(-%.2fz)', max_val, x))
    title('Acoustic wave')
    xlabel('time(ns)')
    ylabel('Pressure')
    ylim([-0.3 2.1])
    saveFig(fig, [fig_path 'p1c.pdf'])
end

%
% Problem 1-d
%

new_ua = 10:10:180;
peak_val = zeros(1, length(new_ua));
for i = 1:length(new_ua)
    time = 0:300;
    pressure = new_ua(i) * GAMMA * H0 * exp(-new_ua(i) * 10^2 * abs(time) * speed * 10^-6);
    [peak, peak_idx] = max(pressure);
    gauss_std = (peak) * 5 / 100;
    pressure = pressure + gauss_std * randn(1, length(pressure));
    [peak, peak_idx] = max(pressure);
    peak_val(i) = peak;
end
func = @(x)(sum(abs(peak_val - x * new_ua)));
x0 = 0;
x = lsqnonlin(func,x0, [], [], options);
if p1d
    fig = figure();
    plot(new_ua, peak_val, 'o', 'linewidth', 2);
    hold on
    plot(new_ua, x * new_ua, 'linewidth', 2);
    legend('peak vs ua', sprintf('y = %fx', x), 'Location', 'NorthWest');
    title('Peak VS ua')
    xlabel('ua');
    ylabel('Peak')
    ylim([-0.1 2.1])
    saveFig(fig, [fig_path 'p1d.pdf'])
end


%
% Problem 1-e
%

new_ua = 10:10:180;
predict_ua = zeros(1, length(new_ua));
for i = 1:length(new_ua)
    time = 0:300;
    pressure = new_ua(i) * GAMMA * H0 * exp(-new_ua(i) * 10^2 * abs(time) * speed * 10^-6);
    [peak, peak_idx] = max(pressure);
    gauss_std = (peak) * 5 / 100;
    pressure = pressure + gauss_std * randn(1, length(pressure));
    
    func = @(x)(sum(abs(pressure - x * GAMMA * H0 * exp(-x * 10^2 * abs(time) * speed * 10^-6))));
    x = lsqnonlin(func, 100, [], [], options);
    predict_ua(i) = x;
end
func = @(x)(sum(abs(predict_ua - x * new_ua)));
x = lsqnonlin(func, 1, [], [], options);
if p1e
    fig = figure();
    plot(new_ua, predict_ua, 'ro', 'linewidth', 3);
    hold on
    plot(new_ua, x * new_ua, 'b--', 'linewidth', 2);
    legend('ua(real) vs ua(curve fit)', sprintf('y = %fx', x), 'Location', 'NorthWest');
    title('ua(real) vs ua(curve fit)')
    xlabel('ua(real)')
    ylabel('ua(curve fit)')
    saveFig(fig, [fig_path 'p1e.pdf'])
end

%
% Problem 1-f
%
fs = 100e6;
fc_lst = [5 10 25 50] * 1e6;
BW_xdc = 0.6;
BWR = -6;
TPE = -40;
impulse = cell(1, length(fc_lst));
for i = 1:length(fc_lst)
    fc_xdc = fc_lst(i);
    tc = gauspuls('cutoff',fc_xdc,BW_xdc,BWR,TPE);
    t  = -tc : 1/fs : tc;
    impulse{i} = gauspuls(t,fc_xdc,BW_xdc);
end

new_ua = 10:10:180;
peak_val = cell(1, length(fc_lst));
predict_ua = cell(1, length(fc_lst));
for i = 1:length(impulse)
    impulse_response = impulse{i};
    peak_val{i} = zeros(1, length(new_ua));
    predict_ua{i} = zeros(1, length(new_ua));
    for j = 1:length(new_ua)
        time = 0:1500;
        pressure = new_ua(j) * GAMMA * H0 * exp(-new_ua(j) * 10^2 * abs(time) * speed * 10^-6);
        [peak, peak_idx] = max(pressure);
        gauss_std = (peak) * 5 / 100;
        pressure = pressure + gauss_std * randn(1, length(pressure));
        %length(pressure)
        pressure = conv(pressure, impulse_response, 'same');
        %close all
        %plot(pressure)
        %dasd
        %length(pressure)
        func = @(x)(sum(abs(pressure - conv(x * GAMMA * H0 * exp(-x * 10^2 * abs(time) * speed * 10^-6), impulse_response, 'same'))));
        predict_ua{i}(j) = lsqnonlin(func, new_ua(j), [], [], options);
        %func(new_ua(j))
        [peak, peak_idx] = max(pressure);
        peak_val{i}(j) = peak;
    end
    %break
end
%plot(pressure)
%sadad
all_x = zeros(1, length(fc_lst));
for i = 1:length(impulse)
    peak_x = peak_val{i};
    func = @(x)(sum(abs(peak_x - x * new_ua)));
    x0 = 180;
    all_x(i) = lsqnonlin(func,x0, [], [], options);
end
if p1f
    fig = figure();
    plot(impulse{1}, 'linewidth', 2);
    xlabel('Frequency')
    ylabel('Amplitude')
    title('Transducer impulse reponse')
    saveFig(fig, [fig_path 'p1f-tx.pdf']);
    
    fig = figure();
    for i = 1:length(impulse)
        subplot(sprintf('22%d', i))
        plot(new_ua, peak_val{i}, 'ro', 'linewidth', 2);
        hold on
        plot(new_ua, all_x(i) * new_ua, '-', 'linewidth', 2);
        legend('peak vs ua', sprintf('y=%fx', all_x(i)), 'Location', 'NorthWest');
        title(sprintf('Peak vs ua ( f = %d MHz)', fc_lst(i) * 10^-6))
        xlabel('ua')
        ylabel('Peak')
    end
    saveFig(fig, [fig_path 'p1f-1.pdf']);
    
    fig = figure();
    for i = 1:length(impulse)
        estimate_ua = predict_ua{i};
        func = @(x)(sum(abs(estimate_ua - x * new_ua)));
        x = lsqnonlin(func, 1, [], [], options);
        
        subplot(sprintf('22%d', i))
        plot(new_ua, estimate_ua, 'ro', 'linewidth', 2);
        hold on
        plot(new_ua, x * new_ua, '-', 'linewidth', 2);
        legend('estimate vs real', sprintf('y=%fx', x), 'Location', 'NorthWest');
        title(sprintf('ua (estimate) vs ua (real) ( f = %d MHz)', fc_lst(i) * 10^-6))
        xlabel('ua (real)')
        ylabel('ua (estimate)')
    end
    saveFig(fig, [fig_path 'p1f-2.pdf']);
end


%
% Problem 2-a and 2-b
%

x = 150;   % in g Hb/liter
MW_Hb = 64500;  % g Hb/mole, molecular weight of Hb
load e_HbO2_Hb; % molar extinction coef. table, (lambda, HbO2, Hb)

SO2 = [0.2 0.4 0.6 0.8 1.0];
lambda = [578 584 590 596];
if p2b
    lambda = [760 780 800 820];
end
ua = zeros(length(SO2), length(lambda));
e_HbO2 = zeros(1, length(lambda));
e_Hb = zeros(1, length(lambda));
peak_val = zeros(length(SO2), length(fc_lst), length(lambda));
predict_SO2 = zeros(length(SO2), length(fc_lst));
for i = 1:length(SO2)
    ua_SO2 = 2.303*e_HbO2_Hb(:,2)*x*SO2(i)/MW_Hb + 2.303*e_HbO2_Hb(:,3)*x*(1-SO2(i))/MW_Hb;
    for j = 1:length(lambda)
        ua(i, j) = ua_SO2(e_HbO2_Hb(:,1) == lambda(j), 1);
        e_HbO2(j) = e_HbO2_Hb(e_HbO2_Hb(:, 1) == lambda(j), 2);
        e_Hb(j) = e_HbO2_Hb(e_HbO2_Hb(:, 1) == lambda(j), 3);
    end
end
%speed_light = 3 * 10^8;
time = 0:1500;
for i = 1:length(SO2)
    for j = 1:length(fc_lst)
        impulse_response = impulse{j};
        for k = 1:length(lambda)
            now_ua = ua(i, k);
            pressure = now_ua * GAMMA * H0 * exp(-now_ua * 10^2 * abs(time) * speed * 10^-6);
            pressure = conv(pressure, impulse_response, 'same');
            [peak, ~] = max(pressure);
            peak_val(i, j, k) = peak;
        end
        A = [e_HbO2' e_Hb'];
        b = reshape(peak_val(i, j, :), 1, 4)';
        x = A \ b;
        predict_SO2(i, j) = x(1) / sum(x);
    end
end

if p2a || p2b
    fig = figure();
    for i = 1:length(fc_lst)
        subplot(sprintf('22%d', i))
        func = @(x)(sum(abs(predict_SO2(:, i)' - x * SO2)));
        slope = lsqnonlin(func, 1.5, [], [], options);
        plot(SO2, predict_SO2(:, i), 'ro', 'linewidth', 2);
        hold on
        plot(SO2, slope * SO2, 'b-', 'linewidth', 2);
        legend('estimated vs real', sprintf('y=%fx', slope), 'Location', 'NorthWest')
        %title('SO2 (estimated) vs SO2 (real)')
        title(sprintf('SO2 (estimate) vs SO2 (real) ( f = %d MHz)', fc_lst(i) * 10^-6))
        xlabel('SO2 (real)')
        ylabel('SO2 (estimated)')
        ylim([0.18 1.1])
    end
    if p2a
        saveFig(fig, [fig_path 'p2a.pdf']);
    elseif p2b
        saveFig(fig, [fig_path 'p2b.pdf']);
    end
end













