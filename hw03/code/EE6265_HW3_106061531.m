%
%
%   EE6265 Fu-En Wang 106061531 HW3 12/10/2017
%
%
clear;

fig_path = '../doc/src/';
GAMMA = 0.12;
H0 = 0.09; % J/cm^2
ua = 180; % cm^-1
speed = 1.5; % mm/us
thick = 300; % um
options = optimset('Display','off');

p1a = false;
p1b = false;
p1d = false;
p1e = true;

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
func = @(x)sum(abs(pressure - max_val*exp(-x * abs(time) * speed * 10^-4)));
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
















