function[time_avg_periodicity,param_avg_periodicity,periods, power_periods] = get_periodicity(path_model, depths_average_for_periodicity, averaged_parameters_for_periodicity,initiate_calculation_peridodicity)

    [time_avg, param_avg] = get_depth_averages(path_model, depths_average_for_periodicity, averaged_parameters_for_periodicity);
    % this is velocity, *100 so to be in cm/yr
    % y = param_avg{:} .* 100;
    for i = 1:numel(param_avg)
        param_avg{i} = param_avg{i} .* 100;
    end
    
    
    for i = 1:numel(param_avg)
        % get data for current parameter
        xxx = time_avg;
        yyy = param_avg{i};
    
        % find periodicity in data after 0.4e9 yrs
        indexes = find(xxx > initiate_calculation_peridodicity./1e6);
%             if isempty(indexes)
%                 % do nothing
%             else
        minimum_time_row = indexes(1);
        maximum_time_row = indexes(end);
        % only take second half of data due to collapse of lid
        xxx = xxx(minimum_time_row:maximum_time_row);
        yyy = yyy(minimum_time_row:maximum_time_row);
%             end
    

    
        % calculate the periodicity of the data
        Fs = 1 / (xxx(2) - xxx(1)); % sampling frequency
        L = length(yyy);
        power_periods = fft(yyy);
        P2 = abs(power_periods / L);
        P1 = P2(1 : L / 2 + 1);
        P1(2:end-1) = 2 * P1(2:end-1);
        f = Fs * (0 : (L/2)) / L;
        periods = 1 ./ f;
    
        f = f(2:end);
        periods = periods(2:end);
    
        % truncate periods and Y to have the same length
        n = min(length(periods), length(power_periods));
        periods = periods(1:n);
        power_periods = power_periods(1:n);
    
    
    end
time_avg_periodicity = time_avg;
param_avg_periodicity   = param_avg;

end