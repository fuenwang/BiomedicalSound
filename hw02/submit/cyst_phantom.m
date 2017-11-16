%
% EE6265 Fu-En Wang 106061531 HW2 11/14/2017
%

function [pos, amp] = cyst_phantom (N, C)

x_size = 15/1000;   %  Width of phantom [mm]
y_size = 0;   %  Transverse width of phantom [mm]
z_size = 20/1000;   %  Height of phantom [mm]
z_start = 30/1000;  %  Start of phantom surface [mm];

%  Creat the general scatterers

lambda = 4.2286e-04; % in m
ggg = 2 * lambda;
grid_x = 0:(ggg):x_size;
grid_y = grid_x;
grid_z = 0:(ggg):z_size;
x = (randsample(grid_x, N, true) / x_size - 0.5)' * x_size;
y = (randsample(grid_y, N, true) - 0.5)' * y_size;
z = (randsample(grid_z, N, true)' / z_size)*z_size + z_start;

%{
rand('state',12345);
x = (rand (N,1)-0.5)*x_size;
y = (rand (N,1)-0.5)*y_size;
z = rand (N,1)*z_size + z_start;
%}
pos=[x y z];
%  Generate the amplitudes with a Gaussian distribution
randn('state',45678);
amp=randn(N,1);

%  Make the cyst and set the amplitudes to zero inside
%  5 mm cyst
r=2.5/1000;      %  Radius of cyst [mm]
%r=2.7/1000;
zc=40/1000;  
xc=0/1000;    %  Place of cyst [mm]

inside = find(((x-xc).^2 + (z-zc).^2) < r^2) ;
amp(inside) = amp(inside)*(10^(C/20));

