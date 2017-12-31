% Template of HW1 Part1 - Baseband Beamformation: Single Element Synthetic Aperture Focusing
%
%							No guarantee
%								Edited by Meng-Lin Li, 12/03/2008
%								Modified by Meng-Lin Li, 12/14/2010
%								Modified by Meng-Lin Li, 12/03/2013
%
%	Phase rotation: to avoid re-modulation, tau for phase rotation = pt2element_delay - 2*r/c where r is the range of the imaging point 
%   tau: relative term over array elements, not absolute term
%								Modified by Meng-Lin Li, 12/20/2013
%								Modified by Meng-Lin Li, Ph. D., 12/16/2015
%								Modified by Meng-Lin Li, Ph. D., 10/27/2016
%								Dept. of Electrical Engineering, National Tsing Hua University

clear all

% ----- (1) Channel Data Simulation ------ %
% --- provide transducer parameters
fc = 5; % transmit center frequency, in MHz
BW = 0.7;	% fractional bandwidth
Nelement = 65; % number of array elements
lambda = soundv/fc; % wavelength, in mm
pitch = lambda/2;  % distance between two adjacent array elements, in mm
                   % you may try  pitch = lambda/4, see if you can suppress the grating lobes

array_ele_coordinate = zeros(Nelement,3);   % cooridnate of each array element, row: element, column: coordinate (x,y,z) in mm
array_ele_coordinate(:,1) = ((1:Nelement) - (Nelement+1)/2)*pitch;


% --- provide target positions
pt_coordinate = [-5 0 10; 0 0 20; 15 0 30];   % coordinate of point scatterers, (x y z) in mm
Npoint = size(pt_coordinate,1);

% --- find the corresponding time delay
soundv = 1.5;  % speed of sound, in mm/us
time_delay = zeros(Nelement, Npoint);
for iPoint = 1:Npoint,
    for iElement = 1:Nelement,
        distance = norm(pt_coordinate(iPoint, :) - array_ele_coordinate(iElement, :));
        time_delay(iElement,iPoint) = 2 * distance / soundv;;        
    end
end    



%  --- emulate analog channel data 
fs_analog = 64*fc; % sampling rate to emulate analog signals, in MHz
tc = gauspuls('cutoff',5E6,.7,-6,-40);
t = -tc : 1/(fs_analog*1E6) : tc;

impulse_response = gauspuls(t,5.0E6,.7); % one way, impulse response of each array element, see HW1 template
impulse_response = impulse_response/max(abs(fft(impulse_response)));
impulse_response_2way = conv(impulse_response,impulse_response); % two way impulse response of each array element

max_time = max(max(time_delay));
Nsample = ceil(max_time*fs_analog);    % number of sample points for a beam or scan-line


channel_data = zeros(Nsample, Nelement); % analog channel data
% locate the delta fuctions
for iElement = 1:Nelement,
    for iPoint = 1:Npoint,
       channel_data(round(time_delay(iElement,iPoint) * fs_analog), iElement) = 1;        % ??? convert time delay to index/ sample points
    end
end    

Npoint_impulse_response_2way = length(impulse_response_2way);
% convolution with impulse_response_2way
for iElement = 1:Nelement,
    tmp = conv(channel_data(:, iElement), impulse_response_2way, 'same'); 
    channel_data(:,iElement) = tmp;
end   

% --- make wavefield plot of analog channel data
???

% --- sampled cahnnel data
fs = ???*fc;	% new sampling rate
D = fs_analog/fs;	% decimation rate, better D is an integer
channel_data = channel_data(???,:); % decimation

% --- make wavefield plot of sampled channel data
???


% ----- (2) Baseband Dynamic Receive beamforming  ---- 
% --- Baseband demodulatoin on channel data, see the similar part of RF beamformer template
% Baseband demodulation: (1) demodulation (2) LPF
% demodulation
[Nsample, Nelement] = size(channel_data);
t = ???;
channel_data = ???; % * exp(-j*2*pi*fc*t), demodulated channel data

% LPF
fcut = ???;
forder = ???; % filter order, determined by checking if the filter response satisfies your filtering requirement
b = fir1(forder,fcut/(fs/2)); % or by fir2, FIR filter design
figure
freqz(b,1); % check filter response

channel_data = conv2(b,1,???,'same'); % baseband channel data

% you may try to further decimate the channel data???
% channel_data = ???;
% fs = ??;
% [Nsample, Nelement] = size(channel_data);


% --- calculate beam spacing and number of beams
dsin_theta = ???; % beam spacing
Nbeam = ???; % number of beams used to sample the 120-degree sector.
w = ???;	% apodization: ones(1,Nelement) or hanning(Nelement)

beam_buffer = zeros(Nsample,Nbeam); % r-sin(theta) beam buffer
% --- Baseband beam formation looping
for iBeam = 1:Nbeam,
    iBeam
    theta = asin(dsin_theta*(iBeam-(Nbeam+1)/2));
    for iSample = 1:NSample,
		r = soundv*(iSample/fs)/2; % the range of the imaging point 
        x_beam = ???;    % x cooridnate of imaging point at (R,sin(theta)), i.e., (iBeam, iSample)
        z_beam = ???;    % z cooridnate of imaging point at (R,sin(theta)), i.e., (iBeam, iSample)
        
        for iElement = 1:Nelement, 
            pt2element_delay = ???; % time delay between imaging point at (R, sin(theta)) or (iBeam, iSample) and Element i
            idx = ???; % convert time delay to "index"
	    phase_rotation = exp(sqrt(-1)*2*pi*fc*(pt2element_delay)); % absolute phase rotation exp(j*2*pi*fd*tau), tau: pt2element_delay, fd: demodulation frequency
		%phase_rotation = exp(sqrt(-1)*2*pi*fc*(pt2element_delay-2*r/soundv)); % relative phase rotation exp(j*2*pi*fd*tau),tau = pt2element_delay - 2*r/c 
            if (idx >= 1) && (idx <= Nsample),
            	% Create beam buffer by computing the coherent sum across the array, or DO SUM FOR EVERY OTHER ELEMENT (33 ELEMENTS SHOULD CONTRIBUTE)
                beam_buffer(iSample, iBeam) = beam_buffer(iSample,iBeam) + w(iElement)*channel_data(idx,iElement)*phase_rotation; % (fine delay and phase rotation) and weighted sum
            end
	
	    % You can create Delayed channel data here ???

        end
        
    end        
end


% --- Display the beam buffer over a logarithmic scale of 40 dB (i.e., 40 dB dynamic range)
DR = 40; % dyanmic range in dB
R_axis = ???
sin_theta_axis = ???
figure
image(R_axis, sin_theta_axis, 20*log10(abs(BBbeam_buffer)/max(max(abs(BBbeam_buffer)))+eps)+DR);
colormap(gray(DR))
colorbar
axis image
xlabel('sin(theta)')
ylabel('R (mm)')
% title(???)

% --- scan conversion
% image
x_size = ?; % can be determined according to the size of scanning area, in mm
z_size = ?;
Npixel_x = ?;	% e.g., 512
Npixel_z = ?;	% e.g., 512
dx = x_size/Npixel_x;
dz = z_size/Npixel_z;

sector_img = zeros(Npixel_z, Npixel_x); % sector image

% scan conversion, with Matlab function interp2() or pcolor()
sector_img = interp2(???, ???, abs(BBbeam_buffer), ???, ???, 'bilinear'); % biliniear interpolation

% display the sector image with 40 dB dynamic range
DR = 40; % dynamic range in dB
x_axis = ???;
z_axis = ???; 
figure
image(x_axis, z_axis, 20*log10(sector_img/max(max(sector_img))+eps)+DR);
colormap(gray(DR))
colorbar
axis image
xlabel('x (mm)')
ylabel('z (mm)')
% title('in dB (envelope detection, DR=40)')


% --- PSF assessment for each point
???




