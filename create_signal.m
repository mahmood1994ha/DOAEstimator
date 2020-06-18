% ## Author: Mahmoud <mahmoud@mahmoud-ThinkPad-E480>
% ## Created: 2020-01-20

function [y t] = create_signal (freq, phase,length,sampling_period)
  if (freq<0 || length<=0 || sampling_period<=0)
    fprintf('error: invalid inputs');
  end
  %create time vector

  period = 1/freq;
  dt = sampling_period;
  fprintf('signal frequency = %f\n',freq);
  disp(period);
  disp(sampling_period);
  %check whether the sampling period is enough 
  if (sampling_period > period/10)
    fprintf('error: not enough samples\n');
    retval = [zeros(length/sampling_period,0) zeros(length/sampling_period,0)];  
  else
    t = (0:dt:length-dt)';
    y = sin(2*pi*freq*t + phase*pi/180);
    retval = [y t];
  end
end
