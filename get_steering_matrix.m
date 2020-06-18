% ## Author: Mahmoud <mahmoud@mahmoud-ThinkPad-X220>
% ## Created: 2020-03-10

function [retval] = get_steering_matrix(sensor_count,signal_angles,sensor_dist,wavelength)
    for i = 1:sensor_count
        A(i,:) = exp(-1i*2*pi*(sensor_dist*(i-1)*sin(signal_angles*pi/180)/wavelength));
    end
    retval = A;
end