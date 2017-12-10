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

p1a = false;
p1b = true;

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
    saveFig(fig, [fig_path 'wave.pdf'])
end

%
% Problem 1-b
%
time = 0:400;
pressure = ua * GAMMA * H0 * exp(-ua * 10^2 * abs(time) * speed * 10^-6);

[max_val, max_idx] = max(pressure);
[min_val, min_idx] = min(pressure);
























