% ## Author: Mahmoud <mahmoud@mahmoud-ThinkPad-X220>
% ## Created: 2020-03-10

function [retval] = bartlett_get_spectrum(R,target_count,theta_samples,sensor_count,sensor_dist,wavelength)
    theta_samples_size = size(theta_samples);
    for theta_scan = 1:theta_samples_size(2)
        %%calculate signal for each radar element
        for sensor = 1:sensor_count
            a(sensor,1) = exp(-1i*2*pi*(sensor_dist*(sensor-1)*sin(theta_samples(theta_scan)*pi/180)/wavelength));
        end
        magnitude = a' * a;
        P_BA(theta_scan,1) = ctranspose(a)*R*a/magnitude; 
    end
    retval = P_BA;
end
