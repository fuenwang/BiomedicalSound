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
save_path = '../doc/src/RF/';
% ----- (1) Channel Data Simulation ------ %
% --- provide transducer parameters
fc = 5; % transmit center frequency, in MHz
BW = 0.7;	% fractional bandwidth
soundv = 1.5;  % speed of sound, in mm/us
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
time_delay = zeros(Nelement, Npoint);
for iPoint = 1:Npoint,
    for iElement = 1:Nelement,
        distance = norm(pt_coordinate(iPoint, :) - array_ele_coordinate(iElement, :));
        time_delay(iElement,iPoint) = 2 * distance / soundv; % in us        
    end
end    


% --- emulate analog channel data 
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
       channel_data(round(time_delay(iElement,iPoint) * fs_analog) ,iElement) = 1;        % ??? convert time delay to index/ sample points
    end
end    

Npoint_impulse_response_2way = length(impulse_response_2way);
% convolution with impulse_response_2way
for iElement = 1:Nelement,
    tmp = conv(channel_data(:, iElement), impulse_response_2way, 'same');   %length of tmp =  length(channel_data(:,iElement))+length(impulse_response_2way) - 1
    %???;   % tmp: from length(channel_data(:,iElement))+length(impulse_response_2way) - 1 => length(channel_data(:,iElement))
    channel_data(:,iElement) = tmp;
end   

% --- make wavefield plot of analog channel data
%???
fig = figure();
imagesc(channel_data);
colorbar;
colormap(gray);
title('channel data (origin)')
saveFig(fig, [save_path '/b-1.pdf'])
%imshow(channel_data)

% --- sampled cahnnel data
fs = 4*fc;	% new sampling rate
D = fs_analog/fs;	% decimation rate, better D is an integer
channel_data = channel_data(1:D:Nsample,:); % decimation

% --- make wavefield plot of sampled channel data
fig = figure();
mx = max(max(channel_data));
imagesc(1:Nelement, 1:D:Nsample, channel_data)
colorbar;
colormap(gray);
title('channel data (sampled)')
saveFig(fig, [save_path '/b-2.pdf'])


%error
% ----- (2) RF Dynamic Receive beamforming  ---- 

% --- calculate beam spacing and number of beams
[Nsample, Nelement] = size(channel_data);


fs_new = 8*fs;
z_axis_old = 0:soundv/(fs):(size(channel_data,1)-1)*soundv/fs;
z_axis_new = 0:soundv/(fs_new):(size(channel_data,1)-1)*soundv/fs;
Nsample_new = length(z_axis_new);
channel_data_new = zeros(length(z_axis_new),size(channel_data,2));
for i = 1 : size(channel_data,2),
    channel_data_new(:,i) = interp1(z_axis_old, channel_data(:,i)', z_axis_new)';
end

dsin_theta = lambda / (2 * (Nelement-1) * pitch); % beam spacing
%dsin_theta = lambda / (1 * (Nelement-1) * pitch); % beam spacing
Nbeam = round(sqrt(3) / dsin_theta); % number of beams used to sample the 120-degree sector.
w = ones(Nelement);	% apodization: ones(1,Nelement) or hanning(Nelement)
%w = hanning(Nelement);
beam_buffer = zeros(Nsample_new,Nbeam); % r-sin(theta) beam buffer

fs = fs_new;
Nsample = Nsample_new;
% --- Note that you need to perform interpolation on acquired channel data here or in the following looping in order to have good enough delay accuracy ---




% --- RF beam formation looping
for iBeam = 1:Nbeam,
    iBeam
    theta = asin(dsin_theta*(iBeam-(Nbeam+1)/2));
    for iSample = 1:Nsample
        z_beam = cos(theta) * (iSample)*soundv/(2*fs_new);
        x_beam = sin(theta) * (iSample)*soundv/(2*fs_new);
        for iElement = 1:Nelement
            %distance = norm(array_ele_coordinate(iElement, :) - [x_beam y_beam z_beam]);
            pt2element_delay = sqrt((abs(x_beam-array_ele_coordinate(iElement,1)).^2+(z_beam-0).^2))/(soundv/2);
            %pt2element_delay = distance / (soundv/2); % time delay between imaging point at (R, sin(theta)) or (iBeam, iSample) and Element i
            idx = round(pt2element_delay * fs); % convert time delay to "index"
            if (idx >= 1) && (idx <= Nsample)
                beam_buffer(iSample, iBeam) = beam_buffer(iSample,iBeam) + w(iElement)*channel_data_new(idx,iElement); % delay and weighted sum
            end
        end
    end        
end
fig = figure();
mx = max(max(beam_buffer));
image(1:Nbeam, 1:Nsample, 255/mx * beam_buffer)
colormap(gray(40))
title('Beam buffer (origin)')
saveFig(fig, [save_path 'b-3-ones.pdf'])

% --- baseband demodulatoin

fig = figure();
[Nsample,Nbeam] = size(beam_buffer);
x_axis = (-fs / 2) : (fs/Nsample) : (fs / 2);
data = beam_buffer(:, round(Nbeam/2));
%data = beam_buffer(:, 30);
plot(x_axis(1:end-1), abs(fftshift(fft(data)))); % find a typical scanline to check the spectrum by fft
xlabel('MHz');
title('FFT of center scanline')
saveFig(fig, [save_path 'b-4.pdf'])

% Baseband demodulation: (1) demodulation (2) LPF
% demodulation
t = ((0:Nsample-1).'/fs)*ones(1,Nbeam);
BBbeam_buffer = beam_buffer.*exp(-(sqrt(-1)*2*pi*fc*t)); % * exp(-j*2*pi*fc*t)

fig = figure(); % check spectrum again
x_axis = (-fs / 2) : (fs/Nsample) : (fs / 2);
data = BBbeam_buffer(:, round(Nbeam/2));
plot(x_axis(1:end-1), abs(fftshift(fft(data))));
xlabel('MHz');
title('FFT of center scanline (demodulation)')
saveFig(fig, [save_path 'b-5.pdf'])

% LPF
fcut = fc;
forder = 80; % filter order, determined by checking if the filter response satisfies your filtering requirement
b = fir1(forder,fcut/(fs/2)); % or by fir2, FIR filter design
fig = figure();
freqz(b,1); % check filter response
title('Frequency response')
saveFig(fig, [save_path 'b-6.pdf'])
%%
BBbeam_buffer = conv2(b,1,BBbeam_buffer,'same'); % baseband data

fig = figure();
mx = max(max(abs(BBbeam_buffer)));
image(1:Nbeam, 1:Nsample, 255/mx * abs(BBbeam_buffer))
colormap(gray(40))
title('Beam buffer (origin2)')

fig = figure(); % check spectrum again
x_axis = (-fs / 2) : (fs/Nsample) : (fs / 2);
data = BBbeam_buffer(:, round(Nbeam/2));
plot(x_axis(1:end-1), abs(fftshift(fft(data))));
xlabel('MHz');
title('FFT of center scanline (demodulation and LPF)')
saveFig(fig, [save_path 'b-7.pdf'])

% --- Display the beam buffer over a logarithmic scale of 40 dB (i.e., 40 dB dynamic range)
DR = 40; % dyanmic range in dB
R_axis = 0:soundv/(2*fs):(Nsample-1)*soundv/(2*fs);
sin_theta_axis = dsin_theta*((1:Nbeam)-(Nbeam+1)/2);
fig = figure();
image(sin_theta_axis, R_axis, 20*log10(abs(BBbeam_buffer)/max(max(abs(BBbeam_buffer)))+eps)+DR);
size(BBbeam_buffer);
colormap(gray(DR))
xlabel('sin(theta)')
ylabel('R (mm)')
title('BBbeam\_buffer (40 dB range)')
saveFig(fig, [save_path 'b-8.pdf'])


% --- scan conversion
% image
x_size = 2*max(sin_theta_axis)*max(R_axis); % can be determined according to the size of scanning area, in mm
z_size = max(R_axis);
Npixel_x = 512;	% e.g., 512
Npixel_z = 512;	% e.g., 512
dx = x_size/(Npixel_x-1);
dz = z_size/(Npixel_z-1);

%sector_img = zeros(Npixel_z, Npixel_x); % sector image
x = (-x_size/2):dx:(x_size/2);
z = 0:dz:z_size;
inter_theta = zeros(Npixel_z, Npixel_x);
inter_R = zeros(Npixel_z, Npixel_x);
for row = 1:Npixel_z
    for col = 1:Npixel_x
      inter_R(row, col) = norm([x(col) z(row)]);
      inter_theta(row, col) = atan(x(col) / z(row));
    end
end
% scan conversion, with Matlab function interp2() or pcolor()
sector_img = interp2(asin(sin_theta_axis) , R_axis, abs(BBbeam_buffer), inter_theta, inter_R, 'linear'); % biliniear interpolation

% display the sector image with 40 dB dynamic range
DR = 40; % dynamic range in dB
x_axis = x;
z_axis = z; 
fig = figure();
image(x_axis, z_axis, 20*log10(sector_img/max(max(sector_img))+eps)+DR);
colormap(gray(DR))
axis image
xlabel('x (mm)')
ylabel('z (mm)')
title('Sector img in dB (envelope detection, DR=40)')
saveFig(fig, [save_path 'b-9.pdf'])

% --- PSF assessment for each point

%}


