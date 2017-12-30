% Template of HW1 Part1 - RF Beamformation: Single Element Synthetic Aperture Focusing
%
%							No guarantee
%								Edited by Meng-Lin Li, 12/03/2008
%								Modified by Meng-Lin Li, 12/14/2010
%								Modified by Meng-Lin Li, 12/03/2013
%								Modified by Meng-Lin Li, Ph. D., 12/16/2015
%								Modified by Meng-Lin Li, Ph. D. 10/27/2016
%								Dept. of Electrical Engineering, National Tsing Hua University

clear all

% ----- (1) Channel Data Simulation ------ %
% --- provide transducer parameters
fc = 5; % transmit center frequency, in MHz
BW = 0.7;	% fractional bandwidth
Nelement = 65; % number of array elements
soundv = 1.5; % mm/us
lambda = soundv/fc; % wavelength, in mm
pitch = lambda/2;  % distance between two adjacent array elements, in mm
                   % you may try  pitch = lambda/4, see if you can suppress the grating lobes


array_ele_coordinate = zeros(Nelement,3);   % cooridnate of each array element, row: element, column: coordinate (x,y,z) in mm
array_ele_coordinate(:,1) = ((1:Nelement) - (Nelement+1)/2)*pitch;


% --- provide target positions
pt_coordinate = [-5 0 10; 0 0 20; 15 0 30];   % coordinate of point scatterers, (x y z) in mm
Npoint = size(pt_coordinate,1);

% --- find the corresponding time delay
%soundv = 1.5;  % speed of sound, in mm/us
time_delay = zeros(Nelement, Npoint);
for iPoint = 1:Npoint,
    zdis = pt_coordinate(iPoint,3);
    x_pos = pt_coordinate(iPoint,1);
    for iElement = 1:Nelement,
        time_delay(iElement,iPoint) = 2*sqrt((array_ele_coordinate(iElement,1)-x_pos)*(array_ele_coordinate(iElement,1)-x_pos)+zdis.^2)/soundv;     
    end
end    



% --- emulate analog channel data 
fs_analog = 64*fc; % sampling rate to emulate analog signals, in MHz

tc = gauspuls('cutoff',5E6,0.7,-6,-40);
t  = -tc : 1/(fs_analog*1E6) : tc;
yi = gauspuls(t,5E6,0.7); plot(t,yi)

impulse_response = yi; % one way, impulse response of each array element, see HW1 template
impulse_response = impulse_response/max(abs(fft(impulse_response)));
impulse_response_2way = conv(impulse_response,impulse_response); % two way impulse response of each array element
%[c,lag] = xcorr(yi,yi);

max_time = max(max(time_delay));
Nsample = ceil(max_time*fs_analog);    % number of sample points for a beam or scan-line (最多需要幾個sample點for 1 scan line)


channel_data = zeros(Nsample, Nelement); % analog channel data
% locate the delta fuctions
for iElement = 1:Nelement,
    for iPoint = 1:Npoint, 
       channel_data(ceil(time_delay(iElement,iPoint)*fs_analog),iElement) = 1; % ??? convert time delay to index/ sample points
    end
end    

Npoint_impulse_response_2way = length(impulse_response_2way);
% convolution with impulse_response_2way
for iElement = 1:Nelement,
    tmp = conv(channel_data(:,iElement), impulse_response_2way, 'same');   %length of tmp = length(tx_sig) + length(channel_data(:,iElement)) - 1, resulted length m+n-1
    %tmp = tmp(length(tx_sig):end);   % tmp: from m+n-1 => length m
    channel_data(:,iElement) = tmp;
end   

% --- make wavefield plot of analog channel data
figure
imagesc(channel_data);
colorbar;
colormap(gray);

% --- sampled cahnnel data
fs = 4*fc;	% new sampling rate
D = fs_analog/fs;	% decimation rate, better D is an integer
channel_data = channel_data(1:D:end,:); % decimation

% --- make wavefield plot of sampled channel data
figure
imagesc(channel_data);
colorbar;
colormap(gray);
%%
% ----- (2) Baseband Dynamic Receive beamforming  ---- 
% --- Baseband demodulatoin on channel data, see the similar part of RF beamformer template
% Baseband demodulation: (1) demodulation (2) LPF
% demodulation
[Nsample, Nelement] = size(channel_data);
t = ((0:Nsample-1).'/fs)*ones(1,Nelement);
channel_data = channel_data.*exp(-(sqrt(-1)*2*pi*fc*t)); % * exp(-j*2*pi*fc*t)

% LPF
fcut = fc;
forder = 100; % filter order, determined by checking if the filter response satisfies your filtering requirement
b = fir1(forder,fcut/(fs/2)); % or by fir2, FIR filter design
figure
freqz(b,1); % check filter response

channel_data = conv2(b,1,channel_data,'same'); % baseband channel data

%%
% you may try to further decimate the channel data???
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% channel_data = ???;
% fs = ??;
% [Nsample, Nelement] = size(channel_data);


% --- calculate beam spacing and number of beams
dsin_theta = lambda/(2*(Nelement-1)*pitch); % beam spacing
Nbeam = ceil(3.^(0.5)/dsin_theta); % number of beams used to sample the 120-degree sector.
w = ones(1,Nelement);	% apodization: ones(1,Nelement) or hanning(Nelement)

beam_buffer = zeros(Nsample,Nbeam); % r-sin(theta) beam buffer
% --- Baseband beam formation looping
for iBeam = 1:Nbeam,
    iBeam
    theta = asin(dsin_theta*(iBeam-(Nbeam+1)/2));
    for iSample = 1:Nsample,
		r = soundv*(iSample/fs)/2; % the range of the imaging point 
        x_beam = r*sin(theta);    % x cooridnate of imaging point at (R,sin(theta)), i.e., (iBeam, iSample)
        z_beam = r*cos(theta);    % z cooridnate of imaging point at (R,sin(theta)), i.e., (iBeam, iSample)
        
        for iElement = 1:Nelement, 
            pt2element_delay = sqrt((abs(x_beam-array_ele_coordinate(iElement,1)).^2+(z_beam-0).^2))/(soundv/2); % time delay between imaging point at (R, sin(theta)) or (iBeam, iSample) and Element i
            idx = ceil(pt2element_delay*fs); % convert time delay to "i
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

%%
% --- Display the beam buffer over a logarithmic scale of 40 dB (i.e., 40 dB dynamic range)
DR = 40; % dyanmic range in dB
R_axis = 0:soundv/(2*fs):(Nsample-1)*soundv/(2*fs);
sin_theta_axis = dsin_theta*((1:Nbeam)-(Nbeam+1)/2);
figure
image(R_axis, sin_theta_axis, 20*log10(abs(beam_buffer)/max(max(abs(beam_buffer)))+eps)+DR);
colormap(gray(DR))
colorbar

xlabel('sin(theta)')
ylabel('R (mm)')
% title(???)
%%
% --- scan conversion
% image
x_size = 2*max(R_axis)*max(sin_theta_axis); % can be determined according to the size of scanning area, in mm
z_size = max(R_axis);
Npixel_x = 512;	% e.g., 512
Npixel_z = 512;	% e.g., 512
dx = x_size/Npixel_x;
dz = z_size/Npixel_z;

x_axis = -x_size/2:dx:x_size/2-dx;
z_axis = 0:dz:511*dz;
[X,Z] = meshgrid(x_axis,z_axis);
%[RR,THETA] = cart2pol(X,Z);
RR = zeros(Npixel_z, Npixel_x);
THETA = zeros(Npixel_z, Npixel_x);
for i = 1:512,
    for j =1:512,
        RR(i,j) = sqrt((abs(X(i,j)).^2+(Z(i,j)).^2));
        THETA(i,j) = atan(X(i,j)/Z(i,j));
    end 
end
THETA(1,257) = 0;%sector_img = zeros(Npixel_z, Npixel_x); % sector image

% scan conversion, with Matlab function interp2() or pcolor()
sector_img = interp2(sin_theta_axis, R_axis, abs(beam_buffer), sin(THETA), RR, 'bilinear');
% display the sector image with 40 dB dynamic range
DR = 40; % dynamic range in dB

figure
image(x_axis, z_axis, 20*log10(sector_img/max(max(sector_img))+eps)+DR);
colormap(gray(DR))
colorbar
axis image
xlabel('x (mm)')
ylabel('z (mm)')
% title('in dB (envelope detection, DR=40)')


% --- PSF assessment for each point





