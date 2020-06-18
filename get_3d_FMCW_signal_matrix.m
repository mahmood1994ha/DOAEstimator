% ## Author: Mahmoud <mahmoud@mahmoud-ThinkPad-X220>
% ## Created: 2020-03-10


%% initial setup
  clc;
  clear;
  close all;
%   pkg load signal;

  %% simulation parameters
  target_count = 3;
  incoming_signal_frequency = 1*1e4;%Hz
  incoming_signal_frontal_phase = [+30 -60 -50];%degrees
  incoming_signal_saggital_phase = [+30 -65 +0];%degrees
  incoming_signal_wavelength = (3*1e8)/incoming_signal_frequency;
  incoming_signal_SNR = 10;
  sensor_dist = incoming_signal_wavelength/2; %satisfies the wavelength condition d<=lambda/2
  sensor_col_count = 10;
  sensor_row_count = 10;
  theta_d = [-90:10:90];
  signal_sample_count = 2*1e4;
  theta_d_size = size(theta_d);
  angle_sample_count = theta_d_size(2);

% test_tau = get_tau_k([30 10],[50 20],[100 50],1,1,1e3,1e3)
%% perform the calculation
X = get_3d_FMCW_signal_matrix_1([1 1 1],incoming_signal_saggital_phase,incoming_signal_frontal_phase,[1 1 1],1e3,1e1,1e-2,4,4,1e-2,1e-2,incoming_signal_SNR);

% X = get_3d_FMCW_signal_matrix_1([1 1], [30 10],[50 20],[1 2],1e3,10,1e-2,4,4,1e-2,1e-2,10);
X_smooth = get_3d_spatially_smoothed_signal(X,2,2,2);
C_out = get_covariance_matrix_from_smoothed_signal(X_smooth);
N = get_noise_subspace_from_covariance_matrix(C_out);

for d_theta = 1:angle_sample_count
    for d_psi = 1:angle_sample_count
        a = get_combined_scanning_steering_vector(theta_d(d_theta),theta_d(d_psi),1,1e6,1e2,1e-3,4,4,sensor_dist,sensor_dist,2,2,2);
        P(d_theta,d_psi) = get_p_music_3d(a,N);
    end
end

mesh(theta_d,theta_d,abs(P));




function [retval] = get_p_music_3d(a,N)
retval = 1/(ctranspose(a) * N * ctranspose(N) * a);
% R_inv = inv(N);
% retval = 1/(ctranspose(a)*R_inv*a); 
end

function [retval] = get_tau_k_p_q(theta,psi,p,q,dx,dy)
c = 3*10^8;
retval = (2/c) * ( (p-1)* dx*cosd(psi)*sind(theta) + (q-1)* dy*sind(psi)*sind(theta));
end 

function [retval] = get_combined_scanning_steering_vector(d_theta,d_psi,d_r,f_c,B,T,M,N,d_x,d_y,l1,l2,l3)
f_s = 2*f_c;
L = 2*f_c*T;
c = 3*10^8;
alpha = B/T;

%create bearing steering vector
for i = 1:M
    for j = 1:N
        a_theta_psi (i,j) = exp(-1i * 2*pi * f_c * get_tau_k_p_q(d_theta,d_psi,i,j,d_x,d_y));
    end
end

a_theta_psi = reshape(a_theta_psi,M*N,1);

%create range steering vector
tau_r = (2/c) * d_r;
for i =  1:L
    a_r (i,1) = exp(1i*2*pi* (f_c*tau_r + alpha*tau_r*i/f_s - 0.5* alpha*tau_r*tau_r));
end

% disp(size(a_theta_psi));
% disp(size(a_r));

a_theta_psi_smooth = a_theta_psi(1:(l1*l2),1);
a_r_smooth = a_r(1:l3,1);

K = kron(a_theta_psi_smooth,a_r_smooth);

retval = K;

end



function [retval] = get_noise_subspace_from_covariance_matrix(C,K)

% [U S N] = svd(C);
% 
% retval = U;
[V,D] = eig(C);
total_eig_count = size(V);
% disp(total_eig_count);
% E_n = V(1:8,1:3);
[E_n sorted_eigv] = sortem(V,D);
E_n = E_n(1:8,4:8);
retval = E_n;
end 



function [retval] = get_covariance_matrix_from_smoothed_signal(X)

% get p1*p2*p3
full_matrix_size = size(X);
p1p2p3  = full_matrix_size(2);
covariance_scaler = 1/(2*p1p2p3);

squared_x = X * ctranspose(X);

J = flip (eye(full_matrix_size(1)));

C = covariance_scaler * (squared_x + J*squared_x'*J);

retval = C;
end

function [retval] = get_3d_spatially_smoothed_signal(X,l1,l2,l3)
%find dimensions of the smoothed matrix
full_matrix_size = size(X);
M = full_matrix_size(1);
N = full_matrix_size(2);
L = full_matrix_size(3);
% find the depth of the smoothed matrix
p1 = M-l1 +1;
p2 = N-l2 +1;
p3 = L-l3 +1;
fprintf('l1 = %d, l2 = %d, l3 = %d \n', l1 , l2 , l3);
fprintf('p1 = %d, p2 = %d, p3 = %d \n', p1 , p2 , p3);

loop_counter = 0;

% while (t<=p3 && m<=p1 && n<=p2)
%    %construct flattned X
%    
% %    fprintf('l1_x = %d, l2_x = %d, l3_x = %d \n', l1 , l2 , l3);
%    
%    X_local = X(m:m+l1 - 1,n:n+l2 - 1,t:t+l3 - 1);
% %    disp(size(X_local));
%    %flatten X_local
%    X_local_flattened = reshape(X_local,l1*l2*l3,1);
%    X_smooth(:,loop_counter) = X_local_flattened;
%    
%    %increment all counters
%    t = t+1;
%    m = m+1;
%    n = n+1;
%    loop_counter = loop_counter+1;
% end

for m = 1 : p1
    for n = 1 : p2
        for t = 1 : p3
            loop_counter  = loop_counter + 1;
            X_local = X(m:m+l1 - 1,n:n+l2 - 1,t:t+l3 - 1);
            X_local_flattened = reshape(X_local,l1*l2*l3,1);
            X_smooth(:,loop_counter) = X_local_flattened;
        end
    end
end

retval = X_smooth;
end





function [retval] = get_3d_FMCW_signal_matrix_1(beta,psi,theta,r,f_c,b,T,M,N,d_y,d_x,snr)
%     sampling frequency is twice the carrier frequency
    f_s = 2*f_c
    L = f_s * T;
    noise_scaler = 1/ (db2mag(snr));
    for p = 1:M
        for q = 1:N
            X(p,q,:) = get_FMCW_scatter_no_noise(beta,psi,theta,r,f_c,b,T,d_y,d_x,p,q) + noise_scaler * randn(1,L);
        end
    end
    retval = X;
end


function [retval] = get_FMCW_scatter_no_noise(beta,psi,theta,r,f_c,b,T,d_y,d_x,p,q)
target_amplitude_size = size(beta);
target_count = target_amplitude_size(2);
f_s = 2*f_c;
L = f_s * T;
alpha = b/T;
l = [1 : 1: L];
x_p_q = zeros(1,L);

for i = 1 : target_count
    tau_k = get_tau_k(psi(i),theta(i),r(i),p,q,d_y,d_x);
    exponent = -1i*2*pi * (f_c * tau_k + (alpha * tau_k / f_s)* l - 0.5 * alpha * tau_k*tau_k);
    x_p_q = x_p_q +  beta(i) * exp(exponent);
end
retval = x_p_q;
end



function tau  = get_tau_k(target_azimuth,target_elevation,target_range,sensor_row_index,sensor_col_index,sensor_row_dist,sensor_col_dist)
if (size(target_azimuth) ~= size(target_elevation) | size(target_range) ~= size(target_elevation) | size(target_range) ~= size(target_azimuth))
    error('azimuth and elevation data not compatible');
end 
c = 3*10^8;
tau = 2./c * (target_range + (sensor_col_index-1).*sensor_col_dist .* cosd(target_azimuth).*sind(target_elevation) + (sensor_row_index - 1) .* sensor_row_dist .* sind(target_azimuth) .* sind(target_elevation));
end



%% this function is used under license from matlab file exchange. See the following link for resource 
% https://www.mathworks.com/matlabcentral/fileexchange/18904-sort-eigenvectors-eigenvalues
function [P2,D2]=sortem(P,D)
% this function takes in two matrices P and D, presumably the output 
% from Matlab's eig function, and then sorts the columns of P to 
% match the sorted columns of D (going from largest to smallest)
% 
% EXAMPLE: 
% 
% D =
%    -90     0     0
%      0   -30     0
%      0     0   -60
% P =
%      1     2     3
%      1     2     3
%      1     2     3
% 
% [P,D]=sortem(P,D)
% P =
%      2     3     1
%      2     3     1
%      2     3     1
% D =
%    -30     0     0
%      0   -60     0
%      0     0   -90


D2=diag(sort(diag(D),'descend')); % make diagonal matrix out of sorted diagonal values of input D
[c, ind]=sort(diag(D),'descend'); % store the indices of which columns the sorted eigenvalues come from
P2=P(:,ind); % arrange the columns in this order
end