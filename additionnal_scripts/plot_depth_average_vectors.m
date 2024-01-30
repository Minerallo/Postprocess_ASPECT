clear all;
clc;close all;
Path ='/Users/ponsm/Desktop/modelblogin/model/globalscale/anulus2d/mobility_function/test/';
depth_for_average = [5000, 15000, 25000]; %m %Should be consistent with one of the depth asked in the aspect prm file
averaged_parameters = {'temperature','cohesions'}; % can be depth temperature adiabatic_temperature adiabatic_pressure adiabatic_density adiabatic_density_derivative velocity_magnitude sinking_velocity rising_velocity friction_angles cohesions yield_stresses viscosity vertical_heat_flux vertical_mass_flux
% Read data from text file
filename = [Path,'depth_average.txt'];
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
if length(depth_for_average) ~= length(intersect(depth_for_average, available_depths))
    error(sprintf('Invalid depth(s) selected. The available depths are:\n%s', num2str(available_depths')));
end


% Plot friction versus time for each depth of 5000 m
for k=1:length(averaged_parameters)
    figure();
depth_colors = lines(length(depth_for_average));
depths = unique(depth);
i=1;
for d = depths'
    if ismember(d, depth_for_average)
        % Indices for current depth
        idx = depth == d;
        % Plot friction versus time for current depth
        plot(time(idx)./1e6, eval([averaged_parameters{k},'(idx)']), 'DisplayName', sprintf('Depth = %d m', d), 'Color', depth_colors(i, :))
        legend_text{i} = sprintf('Depth = %d m', d);
       hold on;
    end
    i=i+1;
end

% Add labels and title to plot
xlabel('Time (My)');
ylabel(strrep(averaged_parameters{k}, '_', ' '));
title(sprintf('Average %s versus Time for Depths of %s m', strrep(averaged_parameters{k}, '_', ' '), num2str(depth_for_average)));

legend(legend_text);
hold off;
end
