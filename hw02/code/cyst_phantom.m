%  Create a computer model of a cyst phantom. 
%  A circular and higher scattering region is embedded in the phantom, 
%  which is located at (x,y,z) = (0, 0, 40) mm and with a diameter of 5 mm
% All scatterers are situated in a box of (x,y,z)=(15,0,20) mm and the box starts 
%  30 mm from the transducer surface.
%
%
%  Calling: [positions, amp] = cyst_phantom (N,C);
%
%  Parameters:  N - Number of scatterers in the phantom
%               C - contrast between the higher scattering region and the background in dB
%  Output:      positions  - Positions of the scatterers.
%               amp        - amplitude of the scatterers.
%
%  Version 1.0, December 19, 1995 by Joergen Arendt Jensen
%                                   
%                                Modified by Meng-Lin Li, 10/08/2008

function [pos, amp] = cyst_phantom (N, C)

x_size = 15/1000;   %  Width of phantom [mm]
y_size = 0;   %  Transverse width of phantom [mm]
z_size = 20/1000;   %  Height of phantom [mm]
z_start = 30/1000;  %  Start of phantom surface [mm];

%  Creat the general scatterers
rand('state',12345);
x = (rand (N,1)-0.5)*x_size;
y = (rand (N,1)-0.5)*y_size;
z = rand (N,1)*z_size + z_start;

pos=[x y z];
%  Generate the amplitudes with a Gaussian distribution
randn('state',45678);
amp=randn(N,1);

%  Make the cyst and set the amplitudes to zero inside
%  5 mm cyst
r=2.5/1000;      %  Radius of cyst [mm]
zc=40/1000;  
xc=0/1000;    %  Place of cyst [mm]

inside = find(((x-xc).^2 + (z-zc).^2) < r^2) ;
amp(inside) = amp(inside)*(10^(C/20));

