  %% initial setup
  clc;
  clear all;
  close all;
%   pkg load signal;

  %% simulation parameters
  target_count = 3;
  incoming_signal_frequency = 1*1e3;%Hz
  incoming_signal_frontal_phase = [+30 -60 -50];%degrees
  incoming_signal_saggital_phase = [+30 -65 +0];%degrees
  incoming_signal_wavelength = (3*1e8)/incoming_signal_frequency;
  incoming_signal_SNR = 20;
  sensor_dist = incoming_signal_wavelength/2; %satisfies the wavelength condition d<=lambda/2
  sensor_col_count = 10;
  sensor_row_count = 10;
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
  n = noise_scaler*randn(sensor_col_count,signal_sample_count);
  
  A_long = get_steering_matrix(sensor_row_count,incoming_signal_frontal_phase,sensor_dist,incoming_signal_wavelength);
  X_long = A_long * y' + n;
  R_long = X_long * X_long' / signal_sample_count;
  P_long = capon_get_spectrum(R_long,target_count,theta_d,sensor_row_count,sensor_dist,incoming_signal_wavelength);

  %latitude DOA
  A_lat = get_steering_matrix(sensor_col_count,incoming_signal_saggital_phase,sensor_dist,incoming_signal_wavelength);
  %create the composite signal coming to each sensor respectively
  X_lat = A_lat * y' + n;

  %calcualte the autocorrelation matrix of the composite signal
  R_lat = X_lat * X_lat';
  R_lat = R_lat/signal_sample_count;

  % calcualte the spatial spectrum of the signal using compared algorithms.
  %P_BA_long = bartlett_get_spectrum(R_long,target_count,theta_d,sensor_col_count,sensor_dist,incoming_signal_wavelength);
  %P_CA_long = capon_get_spectrum(R_long,target_count,theta_d,sensor_col_count,sensor_dist,incoming_signal_wavelength);
  P_lat = capon_get_spectrum(R_lat,target_count,theta_d,sensor_col_count,sensor_dist,incoming_signal_wavelength);
  figure('name','elevation angles');
  plot(theta_d,P_lat,'color','b');
  grid on;
  title('elevation angles');
  axis([-90 90 -inf inf]);

  figure('name','azimuth angles');
  plot(theta_d,P_long,'color','g');
  grid on;
  title('azimuth angles');
  axis([-90 90 -inf inf]);
  
  %combine spectrums
  P_comb = P_lat *P_long' / signal_sample_count ;

  figure;
  mesh(theta_d,theta_d,abs(P_comb));

  %combine the two spatial spectrum by peak finding