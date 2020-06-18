% ## Author: Mahmoud <mahmoud@mahmoud-ThinkPad-X220>
% ## Created: 2020-03-10

function [retval] = MUSIC_get_spectrum(R,target_count,theta_samples,sensor_count,sensor_dist,wavelength)
    [V,D] = eig(R);
    %% find the noise subspace
    E_n = V(:,1:(sensor_count - target_count));
    theta_samples_size = size(theta_samples);
    for theta_scan = 1:theta_samples_size(2)
        %%calculate signal for each radar element
        for sensor = 1:sensor_count
            a(sensor,1) = exp(-1i*2*pi*(sensor_dist*(sensor-1)*sin(theta_samples(theta_scan)*pi/180)/wavelength));
        end
        P_MUSIC(theta_scan,1) = 1/(ctranspose(a)*E_n*ctranspose(E_n)*a); 
    end
    retval = P_MUSIC;
end