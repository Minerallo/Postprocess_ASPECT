function [time_avg,param_avg] = get_depth_averages_origin(path_model, depths_for_average, averaged_parameters)
% Read data from text file
filename = [path_model,'depth_average.txt'];
fileID = fopen(filename,'r');
header = fgetl(fileID); % Read the header line
header = strrep(header, '#', ''); % Remove '#' character from header
allowed_parameters = regexp(header, '\s+', 'split'); % Split header line using one or more spaces as delimiter
fclose(fileID);

% Check if averaged_parameter is allowed
if ~ismember(averaged_parameters, allowed_parameters)
    error(sprintf('Invalid parameter selected. The allowed parameters are:\n%s', strjoin(allowed_parameters, ', ')));
end

% Read data from text file
delimiter = ' ';
formatSpec = '%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%[^\n\r]';
fileID = fopen(filename,'r');
dataArray = textscan(fileID, formatSpec, 'Delimiter', delimiter, 'MultipleDelimsAsOne', true, 'TextType', 'string', 'HeaderLines' ,1, 'ReturnOnError', false);
fclose(fileID);

% Extract data columns
time = dataArray{1};
depth = dataArray{2};
temperature = dataArray{3};
adiabatic_temperature = dataArray{4};
adiabatic_pressure = dataArray{5};
adiabatic_density = dataArray{6};
adiabatic_density_derivative = dataArray{7};
velocity_magnitude = dataArray{8};
sinking_velocity = dataArray{9};
rising_velocity = dataArray{10};
friction_angles = dataArray{11};
cohesions = dataArray{12};
yield_stresses= dataArray{13};
viscosity = dataArray{14};
vertical_heat_flux = dataArray{15};
vertical_mass_flux = dataArray{16};

% Check if depth_for_average contains only depths that are in depth
available_depths = unique(depth);
if length(depths_for_average) ~= length(intersect(depths_for_average, available_depths))
    error(sprintf('Invalid depth(s) selected. The available depths are:\n%s', num2str(available_depths')));
end


% Plot friction versus time for each depth of 5000 m
time_avg = [];
param_avg = cell(1, length(depths_for_average));
depths = unique(depth);
i = 1;
for d = depths'
    if ismember(d, depths_for_average)
        % Indices for current depth
        idx = depth == d;
        % Plot friction versus time for current depth
        time_avg = time(idx)./1e6;
        param_avg{i}= eval([averaged_parameters,'(idx)']);
    end
i = i + 1;
end

end
