%
% EE6265 Fu-En Wang 106061531 HW2 11/14/2017
%

% start up the Field2 simulation system
clear;
addpath('./Field2')
flag = true;
if flag
    field_init(0)

    % --------------------- simulation flag ----------------------------
    flags.phan_type = 1;		% phantom type, 0: wire phantom, 1: tissue mimicking phantom

    % -------- common parameters of excitation(tx pulse) and transducer response --------------
    BWR = -6;			% (in dB), a ref level used to measure BW ,w.r.t. signal peak 
    TPE = -40;			% trailing pulse envelope(in freq) (in dB)
    fs = 100e6;       % sampling frequency (in Hz)
    soundv = 1480;    % sound speed (in m/s)

    % ---------------------- Transducer parameters ----------------------------
    fc_xdc = 3.5e6;               			% transducer center frequency (in Hz)
    BW_xdc = 0.8;									% fractional bandwidth of transducer
    lambda = soundv/fc_xdc;      				% wave length (in m)
    R = 5/1000;             					%  Radius of transducer (in m)
    apersize = 2*R;								% aperture size (in m)
    Rfocus = 40/1000;       					%  Geometric focus point (in m)
    ele_size = lambda/2;   	 					%  Size of mathematical elements (in m)
    fnumber = Rfocus/apersize;					% f number

    % ----------------- excitation(tx pulse) parameters -------------------------
    fc = 3.5e6;							% central frequency of excitation pulse (in Hz)
    BW = 0.33;							% fractional bandwidth of pulse 

    % ---------------- linear scan parameters ----------------------
    dx = lambda/2;	% distance per step (in m)
    Nstep = 151;	% number of steps (odd preferred)

    % --------- set the sampling frequency -----------
    set_field('fs',fs);
    set_field('c', soundv);

    % -------------- generate transducer handler -------------------
    Th = xdc_concave (R, Rfocus, ele_size);

    % -------------- set the impulse response and excitation of the tx and rx aperture -------------
    % generate impulse response of the tx and rx aperture
    tc = gauspuls('cutoff',fc_xdc,BW_xdc,BWR,TPE);
    t  = -tc : 1/fs : tc;
    impulse_response = gauspuls(t,fc_xdc,BW_xdc);

    % generate excitation pulse
    tc = gauspuls('cutoff',fc,BW,BWR,TPE);
    t  = -tc : 1/fs : tc;
    excitation = gauspuls(t,fc,BW);

    xdc_impulse(Th, impulse_response); 	% set the impulse response for the transducer
    xdc_excitation(Th, excitation);		% set the excitation pulse for the transducer


    % ------------------- do linear scan imaging for single crystal element ---------------------------

    % set phantom type
    if (flags.phan_type == 0),
       scatter_pos_ini = [ 0 0 20/1000; 0 0 40/1000; 0 0 60/1000];
       scatter_amp = ones(3,1);
    else
        % 20 ---- 33
        % 10 ---- 14
        %  5 ---- 3.92
       C = 20;   % Contrast: in dB
       [scatter_pos_ini, scatter_amp] = cyst_phantom(5000,C);   % make the phantom
    end
    scatter_pos = scatter_pos_ini;

    Nscatter = size(scatter_pos_ini,1);
    tstart = zeros(1, Nstep);
    for i = 1:Nstep,
        i
        % reset scatter_pos for each step
        scatter_pos(:,1) = scatter_pos_ini(:,1) + (i-(Nstep+1)/2)*dx;

        % calculate the rx signal
       [tmp, tstart(i)] = calc_scat(Th, Th, scatter_pos, scatter_amp);	% received signal
        rf_tmp(1:length(tmp),i) = tmp;

    end

    % align with the min_tstart 
    min_tstart = min(tstart);
    min_sample = min(min_tstart)*fs;
    for i = 1:Nstep,
        tmp = [zeros(round(tstart(i)*fs-min_sample),1);rf_tmp(:,i)];
        rf_data(1:length(tmp),i) = tmp;
    end
    tstart = min_tstart;    % time offset
else
    load data.mat
end


envelope = abs(hilbert(rf_data));   % !!!! envelope, i.e., the amplitude, find the amplitude, E
intensity = envelope.^2;            % !!!! intensity

envelope_lighter = envelope(1350:1800,65:85);
envelope_darker = envelope(320:954,50:100);
intensity_lighter = intensity(1350:1800,65:85);
intensity_darker = intensity(320:954,50:100);

mean_amp_lighter = mean(envelope_lighter(:));
std_amp_lighter = std(envelope_lighter(:));
ratio_amp_lighter = mean_amp_lighter / std_amp_lighter;

mean_intensity_lighter = mean(intensity_lighter(:));
std_intensity_lighter = std(intensity_lighter(:));
ratio_intensity_lighter = mean_intensity_lighter / std_intensity_lighter;

mean_amp_darker = mean(envelope_darker(:));
std_amp_darker = std(envelope_darker(:));
ratio_amp_darker = mean_amp_darker / std_amp_darker;

mean_intensity_darker = mean(intensity_darker(:));
std_intensity_darker = std(intensity_darker(:));
ratio_intensity_darker = mean_intensity_darker / std_intensity_darker;


envelope_lighter_dB = 20*log10(envelope_lighter);
envelope_darker_dB = 20*log10(envelope_darker);
intensity_lighter_dB = 20*log10(intensity_lighter);
intensity_darker_dB = 20*log10(intensity_darker);

%contrast = abs(mean(intensity_lighter_dB(:)) - mean(intensity_darker_dB(:)));
contrast = abs(mean(envelope_lighter_dB(:)) - mean(envelope_darker_dB(:)));
contrast

envelope_dB = 20*log10(envelope/max(max(envelope)));
srcPath = '../doc/src';
%name = sprintf('%s/radius_%d.pdf', srcPath, C);
name = sprintf('%s/wave_10.pdf', srcPath);
%%{
fig = figure();
image(((1:Nstep)-(Nstep+1)/2)*dx*1000, (tstart+(0:size(rf_data,1))/fs)*soundv/2*1000, envelope_dB+40);
%image(1:Nstep, 0:size(rf_data,1), envelope_dB+40);
colormap(gray(40))
axis image
xlabel('Lateral position (mm)')
ylabel('Depth (mm)') 
%saveFig(fig, name);
if flag  
    % free space for apertures
    xdc_free(Th);
    field_end
end
%%}


   


