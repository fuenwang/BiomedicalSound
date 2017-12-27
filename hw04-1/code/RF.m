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
        time_delay(iElement,iPoint) = distance / soundv; % in us        
    end
end    


% --- emulate analog channel data 
fs_analog = 64*fc; % sampling rate to emulate analog signals, in MHz
tc = gauspuls('cutoff',5.0E3,.6,[],-40);
t = -tc : 1E-6 : tc;
impulse_response = gauspuls(t,5.0E3,.6); % one way, impulse response of each array element, see HW1 template
impulse_response_2way = impulse_response.^2; % two way impulse response of each array element

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
%fig = figure();
mx = max(max(channel_data));
%image(1:Nelement, 1:Nsample, 255 / mx * channel_data)
%imshow(channel_data)

% --- sampled cahnnel data
fs = 4*fc;	% new sampling rate
D = fs_analog/fs;	% decimation rate, better D is an integer
channel_data = channel_data(1:D:Nsample,:); % decimation

% --- make wavefield plot of sampled channel data
%fig = figure();
mx = max(max(channel_data));
%image(1:Nelement, 1:D:Nsample, 255 / mx * channel_data)


%error
% ----- (2) RF Dynamic Receive beamforming  ---- 

% --- calculate beam spacing and number of beams
[Nsample, Nelement] = size(channel_data);
dsin_theta = lambda / (2 * (Nelement-1) * pitch); % beam spacing
Nbeam = round(sqrt(3) / dsin_theta); % number of beams used to sample the 120-degree sector.
w = hanning(Nelement);	% apodization: ones(1,Nelement) or hanning(Nelement)
beam_buffer = zeros(Nsample,Nbeam); % r-sin(theta) beam buffer

% --- Note that you need to perform interpolation on acquired channel data here or in the following looping in order to have good enough delay accuracy ---


% --- RF beam formation looping
tmp = zeros(Nelement, 3);
%w_mat = repmat(w',[Nelement 1]);
for iBeam = 1:Nbeam,
    %iBeam
    theta = asin(dsin_theta*(iBeam-(Nbeam+1)/2));
    %rotation = [cos(theta) -sin(theta); sin(theta) cos(theta)];
    loc_z = cos(theta) * (1:Nsample); 
    loc_x = sin(theta) * (1:Nsample);
    loc_y = zeros(1, Nsample);
    for iSample = 1:Nsample,
        %x_beam = ;    % x cooridnate of imaging point at (R,sin(theta)), i.e., (iBeam, iSample)
        %z_beam = ;    % z cooridnate of imaging point at (R,sin(theta)), i.e., (iBeam, iSample)
        %loc = (rotation * [iSample 0]');
        %z_beam = loc(1);
        %x_beam = loc(2);
        z_beam = loc_z(iSample);
        x_beam = loc_x(iSample);
        y_beam = loc_y(iSample);
        
        tmp(:, 1) = (array_ele_coordinate(:, 1) - x_beam).^2;
        tmp(:, 2) = (array_ele_coordinate(:, 2) - y_beam).^2;
        tmp(:, 3) = (array_ele_coordinate(:, 3) - z_beam).^2;
        idx = round(sqrt(sum(tmp')) / soundv * fs);
        cond = idx >= 1 & idx <= Nsample;
        %idx = idx(idx >= 1 & idx <= Nsample);
        idx_flat = Nsample * (0:(Nelement-1)) + idx;
        idx_flat = idx_flat(cond);
        %size(channel_data(idx, :) * w)
        beam_buffer(iSample, iBeam) = sum(channel_data(idx_flat) * w(cond, :));
        %gg
        %{
        for iElement = 1:Nelement, 
            %pt2element_delay = ???; % time delay between imaging point at (R, sin(theta)) or (iBeam, iSample) and Element i
            pt2element_delay = norm(array_ele_coordinate(iElement, :) - [x_beam y_beam z_beam]) / soundv;
            idx = round(pt2element_delay * fs); % convert time delay to "index"

            if (idx >= 1) && (idx <= Nsample),
            	% Create beam buffer by computing the coherent sum across the array, or DO SUM FOR EVERY OTHER ELEMENT (33 ELEMENTS SHOULD CONTRIBUTE)
                beam_buffer(iSample, iBeam) = beam_buffer(iSample,iBeam) + w(iElement)*channel_data(idx,iElement); % delay and weighted sum
                %beam_buffer(iSample, iBeam) = beam_buffer(iSample,iBeam) + 1*channel_data(idx,iElement);
                % You can create Delayed channel data here ???
                
            end

        end
        %}
    end        
end
mx = max(max(beam_buffer));
image(255 / mx * beam_buffer)
%{
% --- baseband demodulatoin

figure
plot(???, ???); % find a typical scanline to check the spectrum by fft
xlabel('MHz');

% Baseband demodulation: (1) demodulation (2) LPF
% demodulation
t = ((0:Nsample-1).'/fs)*ones(1,Nbeams);
BBbeam_buffer = beam_buffer.*exp(-(sqrt(-1)*2*pi*fc*t)); % * exp(-j*2*pi*fc*t)

figure % check spectrum again
plot(???, ???); % by fft
xlabel('MHz');

% LPF
fcut = ???;
forder = ???; % filter order, determined by checking if the filter response satisfies your filtering requirement
b = fir1(forder,fcut/(fs/2)); % or by fir2, FIR filter design
figure
freqz(b,1); % check filter response

BBbeam_buffer = conv2(b,1,BBbeam_buffer,'same'); % baseband data

figure % check spectrum again
plot(???, ???); % fft
xlabel('MHz');

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
sector_img = interp2(?, ?, abs(BBbeam_buffer), ?, ?, 'bilinear'); % biliniear interpolation

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

%}


