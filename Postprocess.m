% This Postprocess file allows to set parameters for the modules as well as to choose which module to call and plot the results.
clear all; close all; clc;

%Crameri colors
addpath './crameri/';

%% Parameter section
% % path_model = '/Path/of/the/model/folder/'

% Example
% path_model = '/Users/ponsm/Desktop/modelblogin/model/globalscale/anulus2d/continents/058_f03_n1e20_UPM_cont_HK2004_umin1e20_drucker/'

%% Path existence check
% Check if path_model ends with "/"
assert(path_model(end) == '/', 'path_model must end with "/"');

if ~exist(path_model, 'dir')
    error(['The directory ' path_model ' does not exist.']);
end

%% Output options
% 1-Output Average parameters depth
% 2-Output statistics
% 3-Output friction angle with Mobility function
% 4-Output RMS sinking velocity for different depths and Periodicity
% 5-Output topography evolution
% 6-Output layer evolution
% 7-Output Dip topography or dip layer evolution

% Model parameters 
%Only required if topo layer is used
model_length =8700e3; % the length is only required if the model is a box
model_height =3000e3;

% Average depth module
% Depths should be consistent with one of the depth asked in the aspect prm file.
% Parameters can be depth temperature adiabatic_temperature adiabatic_pressure adiabatic_density ...
% adiabatic_density_derivative velocity_magnitude sinking_velocity rising_velocity friction_angles ...
% cohesions yield_stresses viscosity vertical_heat_flux vertical_mass_flux

% Average depth module
depths_for_average = [5000, 15000, 25000]; %meters
averaged_parameters = {'temperature','yield_stresses', 'friction_angles'};


% Specify the statistic parameters to plot
%if 'all' then output everything
statistic_parameters = {'RMS','temperature'};
%     {'Mobility', 'RMS','Divergence','Radial RMS','tangential RMS', 'total RMS velocity'};
depth_average_for_friction = 5000;

% Generate a plot combining averaged parameter 'friction angles' and statistic parameters 'Mobility'
Output_combined_mobility_friction = 'true';


% Average a parameter at different depth and calculate the periodicity
% the module is based on the average_depth module
initiate_calculation_peridodicity = 400e6; %Myr
depths_average_for_periodicity = [5000, 15000, 25000, 35000, 45000, 55000, 65000, 75000, 85000, 95000]; %meters
averaged_parameters_for_periodicity = 'sinking_velocity';

% Topography, Dynamic topography and layer topography evolution can be plotted
% in an interval of time in year defined by the user
plot_topography_start = 0;
plot_topography_end = 300e6; %to be changed by the user

% Topography
postprocess_topography = 'false';
dt_topography = 5e6; %should be consistent with the prm file
resample_topography = 1; %Take topography every x files

% Dynamic topography
postprocess_dynamic_topography = 'true';
resample_dynamic_topography = 200;

% Layer topography
%writing everytime step will need to be given an interval
postprocess_topography_layer = 'true';
dt_topography_layer = 1e6; %should be consistent with the prm file
resample_topography_layer = 10; %Take topography layer data every x files
%By default layer elevation is read from top to bottom such as to track the top of a slab subduction
%But one may want the elevation or thickness of a deeper layer such as
%llsvps material in which case the elevation will need to be read from
%bottom to top.
read_layer_elevation_from_bottom_to_top = 'true'; 


%Dip calculation
% calculate the dip of any topography or topolayer
% will only work if postprocess_topography or postprocess_dynamic_topography
% or postprocess_topography_layer are set to true
calculate_topography_dip = 'true';
calculate_topography_layer_dip = 'true';
% Smoothing Interval for Topography
% This parameter controls the spacing between points used for topography smoothing,
% helping reduce noise. The default x-resolution is set to 10,000 points.
% Users should adjust this parameter based on the final time step analysis
% of the topography from the dip smoothing figures to achieve the desired smoothing effect.
topography_smoothing_interval_for_dip_calculation = 100;


% Heat flux
%writing everytime step will need to be given an interval
postprocess_heatflux= 'false';
dt_heatflux = 1e6; %should be consistent with the prm file
resample_heatflux = 5; %Take heatflux every x files




%%% Plot section

% Get the data and corresponding statistic numbers and indices
[data_stats, stats_number, header, stat_indices] = get_statistics(path_model, statistic_parameters);

% Statistics
time = data_stats(:, 2);
if strcmp(statistic_parameters, 'all')
    for i=1:stats_number
        figure();
        plot(time,data_stats(:,i));xlabel('Time [yr]'); ylabel(header(i));title([header(i), 'versus Time']);  %Ma
    end
else
    for i = 1:length(stat_indices)
        figure();
        plot(time, data_stats(:, stat_indices(i)));
        xlabel('Time [yr]');
        ylabel(header(stat_indices(i)));
        title([header(stat_indices(i)), ' versus Time']);
    end
end

try
% Mobility with friction
if strcmp(Output_combined_mobility_friction, 'true')
    % Get the data and corresponding statistic numbers and indices
    [data_stats, stats_number, header, stat_indices] = get_statistics(path_model, {'Mobility'});
    averaged_parameters_save = averaged_parameters;
    depths_for_average_save = depths_for_average;
    depths_for_average = depth_average_for_friction;
    averaged_parameters = {'friction_angles'};
    % Call the get_depth_averages function to get the average friction angles
    [time_avg, param_avg] = get_depth_averages(path_model, depths_for_average, averaged_parameters);

    % Plot the corresponding statistics and friction angles
    figure();
    ax1 = gca;
    plot(time_avg, param_avg{1}, 'DisplayName', 'Friction angles', 'Color', 'r');
    xlabel('Time [yr]');
    ylabel('Friction angles');
    ax2 = axes('Position',get(ax1,'Position'),...
        'YAxisLocation','right',...
        'Color','none');
    hold(ax2, 'on');
    for i = 1:length(stat_indices)
        plot(ax2, data_stats(:, 2), data_stats(:, stat_indices(i)), 'DisplayName', header(stat_indices(i)));
    end
    hold(ax2, 'off');
    ylabel(ax2, 'Average Mobility');
    title('Average Friction angles and Mobility versus Time');
    legend('Location', 'best');
    averaged_parameters = averaged_parameters_save;
    depths_for_average = depths_for_average_save;
end
catch
    disp('No mobility statistics found.');
end


% Average parameter module
try
depth_colors = lines(length(depths_for_average));
for k = 1:length(averaged_parameters)
    % Call the get_depth_averages function to get the average values
    [time_avg, param_avg] = get_depth_averages(path_model, depths_for_average, averaged_parameters{k});
    figure();
    hold on;
    for i = 1:length(depths_for_average)
        plot(time_avg, param_avg{i}, 'DisplayName', sprintf('Depth = %d m', depths_for_average(i)), 'Color', depth_colors(i, :));
    end
    hold off;
    legend('Location', 'best');
    xlabel('Time (My)');
    ylabel(strrep(averaged_parameters{k}, '_', ' '));
    title(sprintf('Average %s versus Time for Depths of %s m', strrep(averaged_parameters{k}, '_', ' '), num2str(depths_for_average)));
end
catch
     disp('No average parameters statistics found.');
end

% Periodicity
try
if(time(end) > initiate_calculation_peridodicity)
    [time_avg_periodicity,param_avg_periodicity,periods, power_periods] = get_periodicity(path_model, depths_average_for_periodicity, averaged_parameters_for_periodicity,initiate_calculation_peridodicity);

    % Plotting parameters
    line_color = {'r', 'g', 'b', 'k', 'm', 'c', 'y', [0.5 0.5 0.5], [0.5 0.5 0]};
    line_style = {'-', '--', '-.', ':'};
    marker = {'o', '+', '*', '.', 'x', 's', 'd', '^', 'v', '>', '<', 'p', 'h'};
    % fig_width = 800;
    % fig_height = 500;

    % Loop through each parameter and plot against time
    % figure('Position', [100, 100, fig_width, fig_height]);
    figure();
    hold on;
    for i = 1:numel(param_avg)
        plot(time_avg_periodicity, param_avg_periodicity{i}, ...
            'Color', line_color{mod(i-1, numel(line_color))+1}, ...
            'LineStyle', line_style{mod(i-1, numel(line_style))+1}, ...
            'Marker', marker{mod(i-1, numel(marker))+1},...
            'LineWidth', 1.5);
    end
    hold off;

    % Set axes labels and title
    xlabel('Time (years)');
    ylabel('Sinking velocity (cm/yr)');
    title('Sinking velocity as a function of time');
    legend(cellfun(@num2str, num2cell(depths_for_average./1000), 'UniformOutput', false));

    % plot periodicity for each param_avg
    figure;
    hold on;
    cmap = jet(numel(param_avg)); % use parula colormap with the same number of colors as the number of depths
    for i = 1:numel(param_avg)
        scatter(periods, abs(power_periods)', 100, 'o', 'filled', 'MarkerEdgeColor', 'k', 'MarkerFaceColor', cmap(i,:));
        xlabel('Period (years)');
        ylabel('Power (x10^3)');
        title(['Periodicity vs Time (Period = ', num2str(periods(2)), ' years)']);
        xlim([0, 600]);
        %     ylim([0, 0.08]);
    end
    legend(cellfun(@num2str, num2cell(depths_for_average./1000), 'UniformOutput', false));
else
    %do not go throught th periodicity
end
catch
    disp('No periodicity statistics found.');
end

%Topography evolution
if strcmp(postprocess_topography, 'true')
try
    [time_elevation_topography, elevation_topography, x_axis_interp, dip_topography] = get_topography_annulus(path_model, dt_topography, resample_topography, calculate_topography_dip, topography_smoothing_interval_for_dip_calculation);
    
    % Find the indices corresponding to the selected time interval
    start_index = find(time_elevation_topography >= plot_topography_start, 1);
    end_index = find(time_elevation_topography <= plot_topography_end, 1, 'last');
    time_elevation_topography=time_elevation_topography(start_index:end_index);
    elevation_topography=elevation_topography((start_index:end_index),:);    
    % Plot the sorted data with adjusted x axis
    figure();
    surf(x_axis_interp,time_elevation_topography,elevation_topography./1e3);shading interp;c=colorbar;demcmap('inc',[5 -8],0.1);ylabel('Time[My]'),xlabel('Annulus Degrees [deg]');zlabel('Elevation[km]');set(gcf,'color','w');
    c.Label.String= "Elevations [km]";set(gcf,'color','w');view(2);%set(gca, 'color', 'none');grid off;set(gca,'XColor', 'none','YColor','none','ZColor','none'); % FaceLighting = 'gour
catch
    disp('Topography files not found. Skipping topography.');
end
end

%Dynamic topography evolution
if strcmp(postprocess_dynamic_topography, 'true')
try
    [time_elevation, elevation, x_axis_interp, dip_topography] = get_dynamic_topography_annulus(path_model, resample_dynamic_topography, calculate_topography_dip, topography_smoothing_interval_for_dip_calculation, time);
    % Find the indices corresponding to the selected time interval
    start_index = find(time_elevation >= plot_topography_start, 1);
    end_index = find(time_elevation <= plot_topography_end, 1, 'last');
    time_elevation=time_elevation(start_index:end_index);
    elevation=elevation((start_index:end_index),:);    
    % Plot the sorted data with adjusted x axis
%     For now let's use the time from statistic but this will have to be change for a resampling time, time_elevation is obsolete.
    figure();
    surf(x_axis_interp,time_elevation,elevation./1e3);shading interp;c=colorbar;demcmap('inc',[5 -8],0.1);ylabel('Time[My]'),xlabel('Annulus Degrees [deg]');zlabel('Elevation[km]');set(gcf,'color','w');
    c.Label.String= "Elevations [km]";set(gcf,'color','w');view(2);%set(gca, 'color', 'none');grid off;set(gca,'XColor', 'none','YColor','none','ZColor','none'); % FaceLighting = 'gour
catch
    disp('Dynamic topography files not found. Skipping dynamic topography.');
end
end

%Layer topography evolution
if strcmp(postprocess_topography_layer, 'true')
try
[time_elevation_layer, elevation_layer, x_axis_interp_layer, dip_layer] = get_topography_layer(path_model, dt_topography_layer, resample_topography_layer, model_length, model_height, calculate_topography_layer_dip, topography_smoothing_interval_for_dip_calculation, read_layer_elevation_from_bottom_to_top);
    % Find the indices corresponding to the selected time interval
    start_index = find(time_elevation_layer >= plot_topography_start, 1);
    end_index = find(time_elevation_layer <= plot_topography_end, 1, 'last'); 
    time_elevation_layer=time_elevation_layer(start_index:end_index);
    elevation_layer=elevation_layer((start_index:end_index),:);
% Plot the sorted data with adjusted x axis
%     For now let's use the time from statistic but this will have to be change for a resampling time, time_elevation is obsolete.
    figure();
 surf(x_axis_interp_layer,time_elevation_layer,elevation_layer./1e3);shading interp;c=colorbar;ylabel('Time[My]'),xlabel('Model Length [km]');zlabel('Elevation[km]');set(gcf,'color','w');
     c.Label.String= "Elevations Layer [km]";set(gcf,'color','w');view(2);%set(gca, 'color', 'none');grid off;set(gca,'XColor', 'none','YColor','none','ZColor','none'); % FaceLighting = 'gour
 catch
     disp('Topography layer files not found.');
end
end

%Calculated dip of topography or layer evolution
if strcmp(calculate_topography_dip, 'true')
    try
        x_axis_dip = x_axis_interp(2:end);
        figure();
        surf(x_axis_dip, time_elevation, dip_topography);
        shading interp;
        c = colorbar;
        % Find the max value of the dip
        max_value = max(max(abs(dip_topography)));
        % Set color limits based on the maximum dip
        clim([-max_value, max_value]);       
        crameri('vik');
        ylabel('Time[My]'), xlabel('Annulus Degrees [deg]');
        zlabel('Dip Topography');
        set(gcf, 'color', 'w');
        c.Label.String = "Dip Topography";
        view(2);
        disp('Dip topography calculation ');

    catch
        disp('Dip topography calculation not asked or not possible.');
    end
elseif strcmp(calculate_layer_dip, 'true')
    try
        x_axis_dip = x_axis_interp(2:end);
        figure();
        surf(x_axis_dip, time_elevation, dip_layer);
        shading interp;
        c = colorbar;
        % Find the max value of the dip
        max_value = max(max(abs(dip_topography)));
        % Set color limits based on the maximum dip
        clim([-max_value, max_value]);       
        crameri('vik');
        ylabel('Time[My]'), xlabel('Annulus Degrees [deg]');
        zlabel('Dip Layer');
        set(gcf, 'color', 'w');
        c.Label.String = "Dip Layer";
        view(2);
        disp('Dip layer calculation ');
    catch
        disp('Dip layer calculation not not asked or not possible.');
    end
end


% Surface heatflux evolution
if strcmp(postprocess_heatflux, 'true')
try
    [time_heatflux, heatfluxmap, x_axis_interp_heatflux] = get_heatflux_annulus(path_model, dt_heatflux, resample_heatflux);
    % Plot the sorted data with adjusted x axis
    figure();
    surf(x_axis_interp_heatflux,time(1:time_heatflux:size(elevation,1)),heatfluxmap.*1e3);shading interp;c = colorbar;crameri('lajolla',12);ylabel('Time[My]'),xlabel('Annulus Degrees [deg]');zlabel('Surface heat flux[W]');caxis([0 120]);set(gcf,'color','w');
    c.Label.String= "Surface heat flux [mW/m]";set(gcf,'color','w');view(2);%set(gca, 'color', 'none');grid off;set(gca,'XColor', 'none','YColor','none','ZColor','none'); % FaceLighting = 'gour
catch
    disp('Heatflux files not found. Skipping heatflux.');
end
end



% Tile all open figures
tilefigs;
