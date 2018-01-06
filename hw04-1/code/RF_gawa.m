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
yi = gauspuls(t,5E6,0.7); 
%%plot(t,yi)

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
title('Channel data wavefield');
xlabel('Channel');
ylabel('Sample points');
% --- sampled cahnnel data
fs = 4*fc;	% new sampling rate
D = fs_analog/fs;	% decimation rate, better D is an integer
channel_data = channel_data(1:D:end,:); % decimation

% --- make wavefield plot of sampled channel data
figure
imagesc(channel_data);
colorbar;
colormap(gray);
title('Sampled channel data wavefield');
xlabel('Channel');
ylabel('Sample points');
%%
% ----- (2) RF Dynamic Receive beamforming  ---- 

% --- calculate beam spacing and number of beams
[Nsample, Nelement] = size(channel_data);
% --- Note that you need to perform interpolation on acquired channel data here or in the following looping in order to have good enough delay accuracy ---
%%%%%%%%%%%%%%%%%%%%%%%%%%%
fs_new = 8*fs;
z_axis_old = 0:soundv/(fs):(size(channel_data,1)-1)*soundv/fs;
z_axis_new = 0:soundv/(fs_new):(size(channel_data,1)-1)*soundv/fs;
Nsample_new = length(z_axis_new);
channel_data_new = zeros(length(z_axis_new),size(channel_data,2));
for i = 1 : size(channel_data,2),
    channel_data_new(:,i) = interp1(z_axis_old, channel_data(:,i)', z_axis_new)';
end
%%%%%%%%%%%%%%%%%%%%%%%%%%%
dsin_theta = lambda/(2*(Nelement-1)*pitch); % beam spacing
%dsin_theta = lambda/((Nelement-1)*pitch); % beam spacing
Nbeam = ceil(3.^(0.5)/dsin_theta); % number of beams used to sample the 120-degree sector.
w = ones(1,Nelement);	% apodization: ones(1,Nelement) or hanning(Nelement)
%w = hanning(Nelement);
beam_buffer = zeros(Nsample_new,Nbeam); % r-sin(theta) beam buffer

%%
% --- RF beam formation looping
for iBeam = 1:Nbeam,
    iBeam
    theta = asin(dsin_theta*(iBeam-(Nbeam+1)/2));
    theta 
    for iSample = 1:Nsample_new,
        x_beam = (iSample)*soundv/(2*fs_new)*sin(theta);    % x cooridnate of imaging point at (R,sin(theta)), i.e., (iBeam, iSample)
        z_beam = (iSample)*soundv/(2*fs_new)*cos(theta);    % z cooridnate of imaging point at (R,sin(theta)), i.e., (iBeam, iSample)
        for iElement = 1:Nelement, 
            pt2element_delay = sqrt((abs(x_beam-array_ele_coordinate(iElement,1)).^2+(z_beam-0).^2))/(soundv/2); % time delay between imaging point at (R, sin(theta)) or (iBeam, iSample) and Element i
            idx = ceil(pt2element_delay*fs_new); % convert time delay to "index"

            if (idx >= 1) && (idx <= Nsample_new),
            	% Create beam buffer by computing the coherent sum across the array, or DO SUM FOR EVERY OTHER ELEMENT (33 ELEMENTS SHOULD CONTRIBUTE)
                %if (ceil(iElement/2)*2 ~= iElement),
                beam_buffer(iSample, iBeam) = beam_buffer(iSample,iBeam) + w(iElement)*channel_data_new(idx,iElement); % delay and weighted sum                
                %end
            end
        end
    end      
end

%%
%for delayed channel data
[a,b] = size(channel_data_new);
delayed_channel_data = zeros(a,b);
point = 3;
x_targetpoint = pt_coordinate(point,1);
z_targetpoint = pt_coordinate(point,3);
complete = 0;
minidx = 10000;
delay_array = zeros(1,b);
for iElement = 1:Nelement,   
    % You can create Delayed channel data here ??? 
    pt2element_delay = sqrt((abs(x_targetpoint-array_ele_coordinate(iElement,1)).^2+(z_targetpoint-0).^2))/(soundv/2); % time delay between imaging point at (R, sin(theta)) or (iBeam, iSample) and Element i
    idx = ceil(pt2element_delay*fs_new);
    delay_array(1,iElement) = idx;
    %channel_data_new(:,iElement) = [channel_data_new(idx:end,iElement);zeros(idx-1,1)];
end   
delay_array_new = delay_array - ones(1,b) * min(delay_array);
for iElement = 1:Nelement,   
    delayed_channel_data(:,iElement) = [channel_data_new(1+delay_array_new(1,iElement):end,iElement);zeros(delay_array_new(1,iElement),1)];
end
figure
imagesc(delayed_channel_data);
colorbar;
colormap(gray);
title('Sampled channel data wavefield');
xlabel('Channel');
ylabel('Sample points');
%%
% --- baseband demodulatoin
%RFdata = beam_buffer(:,ceil(Nbeam/2));
RFfft = fftshift(abs(fft(beam_buffer)));
f_axis = -fs_new/2:fs_new/length(RFfft):fs_new/2-fs_new/length(RFfft);
figure
plot(f_axis, RFfft); % find a typical scanline to check the spectrum by fft
%xlabel('MHz');
title('Spectra of RF beamforming');
xlabel('Frequency(MHz)');
%ylabel('Sample points');
%%
% Baseband demodulation: (1) demodulation (2) LPF
% demodulation
t = ((0:Nsample_new-1).'/fs_new)*ones(1,Nbeam);
BBbeam_buffer = beam_buffer.*exp(-(sqrt(-1)*2*pi*fc*t)); % * exp(-j*2*pi*fc*t)
BBRFdata = BBbeam_buffer(:,ceil(Nbeam/2));
BBRFfft = fftshift(abs(fft(BBbeam_buffer)));
BBf_axis = -fs_new/2:fs_new/length(BBRFfft):fs_new/2-fs_new/length(BBRFfft);
figure % check spectrum again
plot(BBf_axis, BBRFfft); % by fft
xlabel('MHz');

% LPF
fcut = fc;
forder = 100; % filter order, determined by checking if the filter response satisfies your filtering requirement
b = fir1(forder,fcut/(fs_new/2)); % or by fir2, FIR filter design
figure
freqz(b,1); % check filter response


BBbeam_buffer = conv2(b,1,BBbeam_buffer,'same'); % baseband data
BBRFdata = BBbeam_buffer(:,ceil(Nbeam/2));
BBRFfft = fftshift(abs(fft(BBbeam_buffer)));
BBf_axis = -fs_new/2:fs_new/length(BBRFfft):fs_new/2-fs_new/length(BBRFfft);
figure % check spectrum again
plot(BBf_axis, BBRFfft); % by fft
xlabel('MHz');
%%
% --- Display the beam buffer over a logarithmic scale of 40 dB (i.e., 40 dB dynamic range)
DR = 40; % dyanmic range in dB
R_axis = 0:soundv/(2*fs_new):(Nsample_new-1)*soundv/(2*fs_new);
%R_axis = 0:soundv/fs:110*soundv/fs;
theta_axis = asin(dsin_theta*((1:Nbeam)-(Nbeam+1)/2));
sin_theta_axis = dsin_theta*((1:Nbeam)-(Nbeam+1)/2);
figure
image(sin_theta_axis, R_axis,20*log10(abs(BBbeam_buffer)/max(max(abs(BBbeam_buffer)))+eps)+DR);
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
%%
%[THETA,RR] = meshgrid(theta_axis,R_axis);
%[XX,YY] = pol2cart(THETA,RR);
%[X,Z] = meshgrid(x_axis,z_axis);
% scan conversion, with Matlab function interp2() or pcolor()
%sector_img = interp2(XX,YY, abs(BBbeam_buffer), X,Z, 'bilinear'); % biliniear interpolation

%sector_img = interp2(sin(theta_axis), R_axis, abs(BBbeam_buffer), sin(theta), rho, 'bilinear');
sector_img = interp2(theta_axis, R_axis, abs(BBbeam_buffer), THETA, RR, 'bilinear');
% display the sector image with 40 dB dynamic range
DR = 40; % dynamic range in dB
%%
figure
image(x_axis, z_axis, 20*log10(sector_img/max(max(sector_img))+eps)+DR);
colormap(gray(DR))
colorbar
axis image
xlabel('x (mm)')
ylabel('z (mm)')
title('Sector image in dB (envelope detection, DR=40)')

%%
% --- PSF assessment for each point

projection = max(20*log10(abs(sector_img)/max(max(abs(sector_img)))+eps));
[pks,locs] = findpeaks(projection,'MinPeakHeight',-5);
res6locs = zeros(1,6);
res6pks = [pks(1)-6 pks(1)-6 pks(2)-6 pks(2)-6 pks(3)-6 pks(3)-6];
res20locs = zeros(1,6);
res20pks = [pks(1)-20 pks(1)-20 pks(2)-20 pks(2)-20 pks(3)-20 pks(3)-20];
number1 = 1;
number2 = 1;
for i = 1:length(projection)-1,
    if (number1<7)
        if (projection(i) < res6pks(number1) && projection(i+1) > res6pks(number1))
            if(ceil(number1/2)*2 ~= number1)
                res6locs(number1) = i;
                number1 = number1+1;
            end
        elseif (projection(i) > res6pks(number1) && projection(i+1) < res6pks(number1))
            if(ceil(number1/2)*2 == number1)
                res6locs(number1) = i;
                number1 = number1+1;
            end
        end
    end
    if (number2<7)
        if (projection(i) < res20pks(number2) && projection(i+1) > res20pks(number2))
            if(ceil(number2/2)*2 ~= number2)
                res20locs(number2) = i;
                number2 = number2+1;
            end
        elseif (projection(i) > res20pks(number2) && projection(i+1) < res20pks(number2))
            if(ceil(number2/2)*2 == number2)
                res20locs(number2) = i;
                number2 = number2+1;
            end
        end
    end
end
latres6 = dx*[res6locs(2)-res6locs(1) res6locs(4)-res6locs(3) res6locs(6)-res6locs(5)]
latres20 = dx*[res20locs(2)-res20locs(1) res20locs(4)-res20locs(3) res20locs(6)-res20locs(5)]
figure
plot(x_axis,projection,locs*dx-x_size/2,pks,'o',res6locs*dx-x_size/2,res6pks,'bo', res20locs*dx-x_size/2,res20pks,'go');
grid on
title('Max axial projection in dB')

%%
projection2 = max(20*log10(abs(sector_img')/max(max(abs(sector_img)))+eps));
[pks2,locs2] = findpeaks(projection2,'MinPeakHeight',-5);
res6locs2 = zeros(1,6);
res6pks2 = [pks2(1)-6 pks2(1)-6 pks2(2)-6 pks2(2)-6 pks2(3)-6 pks2(3)-6];
res20locs2 = zeros(1,6);
res20pks2 = [pks2(1)-20 pks2(1)-20 pks2(2)-20 pks2(2)-20 pks2(3)-20 pks2(3)-20];
number3 = 1;
number4 = 1;
for i = 1:length(projection2)-1,
    if (number3<7)
        if (projection2(i) < res6pks2(number3) && projection2(i+1) > res6pks2(number3))
            if(ceil(number3/2)*2 ~= number3)
                res6locs2(number3) = i;
                number3 = number3+1;
            end
        elseif (projection2(i) > res6pks2(number3) && projection2(i+1) < res6pks2(number3))
            if(ceil(number3/2)*2 == number3)
                res6locs2(number3) = i;
                number3 = number3+1;
            end
        end
    end
    if (number4<7)
        if (projection2(i) < res20pks2(number4) && projection2(i+1) > res20pks2(number4))
            if(ceil(number4/2)*2 ~= number4)
                res20locs2(number4) = i;
                number4 = number4+1;
            end
        elseif (projection2(i) > res20pks2(number4) && projection2(i+1) < res20pks2(number4))
            if(ceil(number4/2)*2 == number4)
                res20locs2(number4) = i;
                number4 = number4+1;
            end
        end
    end
end
axires6 = dx*[res6locs2(2)-res6locs2(1) res6locs2(4)-res6locs2(3) res6locs2(6)-res6locs2(5)]
axires20 = dx*[res20locs2(2)-res20locs2(1) res20locs2(4)-res20locs2(3) res20locs2(6)-res20locs2(5)]

figure
plot(z_axis,projection2,locs2*dz,pks2,'o',res6locs2*dz,res6pks2,'bo',res20locs2*dz,res20pks2,'go');
title('Max lateral projection in dB')
grid on


