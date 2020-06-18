  %% initial setup
  clc;
  clear all;
  close all;
%   pkg load signal;

  %% simulation parameters
  target_count = 3;
  incoming_signal_frequency = 1*1e3;
  incoming_signal_phase = [+30 -60 -50];%degrees
  incoming_signal_wavelength = (3*1e8)/incoming_signal_frequency;
  incoming_signal_SNR = 20;
  sensor_dist = incoming_signal_wavelength/2; %satisfies the wavelength condition d<=lambda/2
  sensor_count = 10;
  theta_d = [-90:0.1:90];
  signal_sample_count = 2*1e4;
  theta_d_size = size(theta_d);
  angle_sample_count = theta_d_size(2);

  %% create signals
  %% the signals coming from different targets have to be uncorrelted.
  [y1 t1] = create_signal(0.7*incoming_signal_frequency,90,1,1/(signal_sample_count));
  [y2 t2] = create_signal(0.6*incoming_signal_frequency,90,1,1/(signal_sample_count));
  [y3 t3] = create_signal(0.5*incoming_signal_frequency,90,1,1/(signal_sample_count));
  y = [y1 y2 y3];

  %  create a white noise vector
  noise_scaler = 1/ (db2mag(incoming_signal_SNR));
  n = noise_scaler*randn(sensor_count,signal_sample_count);
  
  %% create a steering matrix (vector for one target)
  for i = 1:sensor_count
    A(i,:) = exp(-1i*2*pi*(sensor_dist*(i-1)*sin(incoming_signal_phase*pi/180)/incoming_signal_wavelength));
  end

  %create the composite signal coming to each sensor respectively
  X = A * y' + n;

  %calcualte the autocorrelation matrix of the composite signal
  R = X * X';
  R = R/signal_sample_count;

  % calcualte the spatial spectrum of the signal using compared algorithms.
  P_BA = bartlett_get_spectrum(R,target_count,theta_d,sensor_count,sensor_dist,incoming_signal_wavelength);
  P_CA = capon_get_spectrum(R,target_count,theta_d,sensor_count,sensor_dist,incoming_signal_wavelength);
  P_MUSIC = MUSIC_get_spectrum(R,target_count,theta_d,sensor_count,sensor_dist,incoming_signal_wavelength);
  
  figure('name','spatial spectrum - comparison between DOA algorithms');

  subplot(3,1,1);
  plot(theta_d,P_BA,'color','r');
  grid on;
  title('Bartlett beamformer');
  axis([-90 90 -inf inf]);

  subplot(3,1,2);
  plot(theta_d,P_CA,'color','b');
  grid on;
  title('CAPON beamformer');
  axis([-90 90 -inf inf]);

  subplot(3,1,3);
  plot(theta_d,P_MUSIC,'color','g');
  grid on;
  title('MUSIC DOA estimation');
  axis([-90 90 -inf inf]);
