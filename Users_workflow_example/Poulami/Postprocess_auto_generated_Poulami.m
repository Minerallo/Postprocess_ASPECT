% This Postprocess file allows to set parameters for the modules as well as to choose which module to call and plot the results.
clear all; close all; clc;

% Paths to add (obsolete)
% addpath(genpath('./Individual_scripts/'))
% addpath(genpath('./Get_scripts/'))
% addpath(genpath('./Import_scripts/'))
% addpath(genpath('./Output_folder/'))
% addpath(genpath('./Plot_scripts/'))

%% Example of output options
% 1-Output Average parameters depth
% 2-Output statistics
% 3-Output friction angle with Mobility function
% 4-Output RMS sinking velocity for different depths and Periodicity
% 5-Output topography evolution
% 6-Output layer evolution
% 7-Output Dip topography or dip layer evolution
% 8-Output depth_average plots
% 9-Output dynamic topography
% 10-Run Paraview spherical postprocess

current_file_path = pwd;
parent_directory = fileparts(current_file_path);
addpath(genpath(parent_directory));
postprocess_path = pwd;

%% Parameter section
% path_model = '/Path/of/the/model/folder/'

% This path should stay one line for the auto postprocess so the line can
% be replaced by a list of models
path_models = {['/Volumes/Jerry/global_models_3d_extract/14_test_250Ma_3D_test_with_shear_heat_01/']};

% Select a folder path where the images and data will be outputted if
% necessary otherwise the default will be './output_folder/'
selected_output_folder_path ='/Users/ponsm/Nextcloud/group_monitoring_earth_evolution_through_time/Research/Michael_Pons/models/Global_model_3D/';
% '/Users/ponsm/Library/CloudStorage/GoogleDrive-mikpons2@gmail.com/My Drive/Postdoc/collab/Rift_copy/Postprocess/';

model_geometry = 'sphere';
% 'box'; 

Save_figures = 'true';
Display_figures = 'false';

% Model parameters

% Average depth module
% Depths should be consistent with one of the depth asked in the aspect prm file.
% Parameters can be depth temperature adiabatic_temperature adiabatic_pressure adiabatic_density ...
% adiabatic_density_derivative velocity_magnitude sinking_velocity rising_velocity friction_angles ...
% cohesions yield_stresses viscosity vertical_heat_flux vertical_mass_flux

% Average depth module
depths_for_average = [5000, 15000, 25000]; %meters
averaged_parameters = {'temperature','yield_stresses', 'friction_angles','viscosity'};
%Average parameter profiles will make some profile evolution of each
%averaged parameters
output_all_averaged_parameters_profiles='true';

% Specify the statistic parameters to plot
%if 'all' then output everything
statistic_parameters = {'all'};
% {'RMS velocity','temperature','Mobility','total RMS velocity'};
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
postprocess_topography = 'true';
dt_topography = 20000; %should be consistent with the prm file
resample_topography = 1; %Take topography every x files
%if a shift should be applied for the topography since the reference frame is not the sea level or geoid but the mantle
% topography_shift = 2400; %m

% Dynamic topography
postprocess_dynamic_topography = 'true';
resample_dynamic_topography = 500;

% Layer topography
%writing everytime step will need to be given an interval
postprocess_topography_layer = 'true';
dt_topography_layer = 1e6; %should be consistent with the prm file
resample_topography_layer = 4; %Take topography layer data every x files
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

%Only required if topo layer is used
model_length =2000e3; % the length is only required if the model is a box
model_height =400e3;

% Heat flux
%writing everytime step will need to be given an interval
postprocess_heatflux= 'false';
dt_heatflux = 1e6; %should be consistent with the prm file
resample_heatflux = 5; %Take heatflux every x files

%% Output additional statistics on geodynamics feature for spherical models
% These statistics require an intermediary step of extraction of the data from paraview using
% the python paraview script Global_extract_series.py accessible in
% Python_paraview
run_spherical_additional_postprocess = 'true';

%Path where the data extracted by paraview are stored locally
path_extracted_files_input ='/Volumes/Jerry/global_models_3d_extract/';

% This process can take a significant amount of time depending on the amount of files to process
%therefore I had a restart option so the postprocessing does not need to
%restart from the first file if it is stopped. 'auto' means that it will
%check first how many steps have been already postprocessed
restart_additional_postprocess = 'auto';

% The user can decide whether or not to output the map projection of the model.
% This can be time-consuming, but the maps can later be used to create time-based videos.
output_figures_for_spherical_additional_postprocess = 'true';

%You probably want to write down the statistics of the additional postprocess
write_geofeatures_statistics = 'true' ;
% You can select which statistics you want to write down
% e.g additional_postprocesses = {'subduction_and_plume_statistics','oceanic_age_statistics','continents_VRMS'};
additional_postprocesses = {'subduction_and_plume_statistics'};

% Select which figures to output
% e.g output_additional_maps_figures = {'topography','geofeatures','strain_rate'};
output_additional_maps_figures = {'geofeatures','strain_rate'};

% In case oceanic age statistics are set, the oceanic age statistics can be calculated without the subduction zones
% This will be slightly more time comsumming as not only the lithosphere data but also the surface data will need to be read too.
remove_subductions_to_oceanic_age = 'true';

% If you use 'subduction_and_plume_statistics' then you may want to tune the tracking of the subduction and plumes
%Select which depths layers number with which you want to process and get some statistics
%from for the plumes [1,..,n]
% For example global_extract series.py, by default extract the parameters
% for the depths 440, 660, 1000, 2000, 2600 km, so the default would be that 1 is shallowest depth layer and 5
% is the deepest

plumes_depths_tracking=[150,440, 1100, 2600, 2800];
%Select the depth layers number for which you want to track subduction,
%alternatively, you can choose 0, so the subductions will be tracked from
%the topography, in that case select at which elevation should a trench be recognized.
subduction_depths_tracking=[150, 440, 1100, 2600, 2800];
plumes_non_adiabatic_tracking_temperature = [50,50,100,250,600];% Non adiabatic threshold for plumes for each depths layer
subduction_non_adiabatic_tracking_temperature = [-100,-100,-100,-100,-100];%% Non adiabatic threshold for subduction for each depths layer
trenches_elevation_threshold = -6000; %elevation threshold for trenches

% In case the user want to plot the gplate distribution of continents
%by default it will plot the distribution or Merdith et al., 2019
plot_continents_border_from_reconstruction='true';

% if the gplate continents distribution is plotted then you need to select
% a geological time at which your model should start.
reference_time_My_to_Ma = 250; 

% The user can load additional field to display as map such as
% different compositions over depth
additional_fields_to_load = {'llsvps','melt_fraction'}; 

additional_fields_threshold = {'>0.05','>0.05','>0.05','>0.25','>0.25';
                                 '','','','',''};

additional_fields_depths_to_visualize = {'150','440','1100','2600','2800';
                                         '150','440','','','' };

%this threshold will only be used to draw the contour of a field selected
%if additional_fields_threshold is set to {''};
default_threshold_fields = 0.25;

% Since the topography in Aspect is relatif to the mantle, we need to make
% it relative to the sea level to compare, therefore we reduce the whole
% topography by 2400 which is the average elevation of the mantle in contact with the surface (ridge)
topography_correction = -2400;

% If the model has continents to track, enter the name of the
% compositional field. By default, continent is the name continent as in the
% Global_extract_series.py and in Prm. If the continents have a different name in the prm, it should also be
% replaced in the Python script.
% Specify which name you use here for the plots.
% Continents_name= 'nan';
compositional_field_name_of_continents = 'continent';

% The Time will be read from the paraview python files that were extracting but in case no Time field was extracted or found
% then we will assign a default output time interval, similarly to the one
% that was given in ASPECT for the vtu output
% Obviously it is much better if the time can be read directly from the
% file as the vtu files outputted by Aspect are not fully fitting the time
% given, it is not a big deal for long term model, but might be ennoying
% for really short one of the order of few tens of thousands of years. 
% In the case where everytime irregualar step is outputted from ASPECT and
% extracted from Paraview then time should be added to the files for consistency. 
Interval_of_time_output_for_additional_postprocess = 1; %in Myr

%% Plot section
% Path existence check
for i=1:numel(path_models)
    % Check if path_model ends with "/"
    assert(strcmp(path_models{i}(end), '/'), 'path_model must end with "/"');

    % Convert cell array to string
    path_model_str = char(path_models{i});

    if ~exist(path_model_str, 'dir')
        error(['The directory ' path_model_str ' does not exist.']);
    end
end
cd(postprocess_path);

% Extracting model titles
model_titles = cell(size(path_models));
for t = 1:numel(path_models)
    path_model = path_models{t};
    %% Creation of the Output repository
    % Split the path into parts
    parts = strsplit(path_model, '/');
    % Extract the last part
    model_title = parts{end-1};

    % Define the output folder
    % Check if selected_output_folder_path is defined
    if exist('selected_output_folder_path', 'var') == 1
        % Use the defined value
        output_folder_path = selected_output_folder_path;
    else
        % Set a default value
        output_folder_path =  './Output_folder/';
    end

    % Create the full path for the new repository
    path_model_output = char(fullfile(output_folder_path, model_title));

    % Check if the repository already exists
    if exist(path_model_output, 'dir') == 7
        fprintf('The repository %s already exists.\n', model_title);
    else
        % Create the new repository
        mkdir(path_model_output);
        fprintf('Repository %s created successfully.\n', model_title);
    end

    %For the additional spherical refinement we want the path where the data layers are stored
    model_title_withbar= [model_title '/'];
    path_model_input = char(fullfile(path_extracted_files_input, model_title_withbar));


    % Get the data and corresponding statistic numbers and indices
    [data_stats, stats_number, header, stat_indices] = get_statistics(path_model, statistic_parameters);


    if strcmp(Save_figures, 'true')
        cleaned_headers = string.empty;  % Initialize an empty string array
        for i=1:numel(header)
            % For saving we clean the header to name the figures properly
            % Remove leading hyphens, consecutive hyphens, special characters, and replace spaces with underscores
            cleaned_header = regexprep(header(i), '^[-]+', '');  % Remove leading hyphens
            cleaned_header = regexprep(cleaned_header, '[-]+', '_');  % Remove consecutive hyphens and replace with a single underscore
            cleaned_header = regexprep(cleaned_header, '[^a-zA-Z0-9_]+', '_');  % Remove special characters
            cleaned_header = regexprep(cleaned_header, '\([^)]*\)', '');  % Remove parentheses and their contents
            cleaned_header = regexprep(cleaned_header, '^_', '');
            cleaned_header = regexprep(cleaned_header, '_+$', '');


            % Replace forward slashes with the word "per"
            cleaned_header = regexprep(cleaned_header, '/', 'per');
            cleaned_headers(end + 1) = cleaned_header;
            % disp(cleaned_header);
        end
    end

    % Statistics
    time = data_stats(:, 2);
    if strcmp(statistic_parameters, 'all')
        for i=1:stats_number
            h=figure;
            plot(time,data_stats(:,i),'r-', 'LineWidth', 2);xlabel('Time [yr]'); ylabel(header(i));title([header(i), 'versus Time']);  %Ma
            set(gcf, 'color', 'w');

            % Set visibility based on Display_figures
            if strcmp(Display_figures, 'true')
                set(h, 'Visible', 'on');
            else
                set(h, 'Visible', 'off');
            end

            if strcmp(Save_figures, 'true')
               save_directory = fullfile(path_model_output, 'statistics_plots');
            if ~exist(save_directory, 'dir')
                mkdir(save_directory);
            end

                fig_filename = fullfile(save_directory, sprintf('Stats_%s.png', cleaned_headers(i)));
                saveas(gcf, fig_filename);
                fprintf('Figure saved: %s\n', fig_filename);
            end

        end
    else
        for i = 1:length(stat_indices)
            h=figure;
            plot(time, data_stats(:, stat_indices(i)),'r-', 'LineWidth', 2);
            set(gcf, 'color', 'w');
            xlabel('Time [yr]');
            ylabel(header(stat_indices(i)));
            title([header(stat_indices(i)), ' versus Time']);

            % Set visibility based on Display_figures
            if strcmp(Display_figures, 'true')
                set(h, 'Visible', 'on');
            else
                set(h, 'Visible', 'off');
            end

            if strcmp(Save_figures, 'true')

            save_directory = fullfile(path_model_output, 'statistics_plots');
            if ~exist(save_directory, 'dir')
                mkdir(save_directory);
            end
    
                fig_filename = fullfile(save_directory, sprintf('Stats_%s.png', cleaned_headers(stat_indices(i))));
                saveas(gcf, fig_filename);
                fprintf('Figure saved: %s\n', fig_filename);
            end
        end
    end

    try
        % Mobility with friction
        if strcmp(Output_combined_mobility_friction, 'true')
            % Get the data and corresponding statistic numbers and indices
            %         [data_stats, stats_number, header, stat_indices] = get_statistics(path_model, {'Mobility'});

            %         %Check index Mobility in statistic_parameters
            index_of_mobility_in_statistics_parameters = find(contains(statistic_parameters, 'Mobility'),1);
            if isempty(index_of_RMS_velocity_in_statistics_parameters)
                disp('Could not plot the mobility with friction as "Mobility" is missing in statistic_parameters')
            end

            averaged_parameters_save = averaged_parameters;
            depths_for_average_save = depths_for_average;
            depths_for_average = depth_average_for_friction;
            averaged_parameters = {'friction_angles'};
            % Call the get_depth_averages function to get the average friction angles
            [time_avg, param_avg] = get_depth_averages(path_model, depths_for_average, averaged_parameters);

            % Plot the corresponding statistics and friction angles
            h=figure;
            ax1 = gca;
            plot(time_avg, param_avg{1}, 'DisplayName', 'Friction angles', 'Color', 'r');
            xlabel('Time [yr]');
            ylabel('Friction angles');
            ax2 = axes('Position',get(ax1,'Position'),...
                'YAxisLocation','right',...
                'Color','none');
            hold(ax2, 'on');
            %         for i = 1:length(stat_indices)
            plot(ax2, time, data_stats(:, stat_indices(index_of_mobility_in_statistics_parameters)), 'DisplayName', header(stat_indices(index_of_mobility_in_statistics_parameters)));
            %         end
            hold(ax2, 'off');
            ylabel(ax2, 'Average Mobility');
            title('Average Friction angles and Mobility versus Time');
            legend('Location', 'best');
            averaged_parameters = averaged_parameters_save;
            depths_for_average = depths_for_average_save;

            % Set visibility based on Display_figures
            if strcmp(Display_figures, 'true')
                set(h, 'Visible', 'on');
            else
                set(h, 'Visible', 'off');
            end

            if strcmp(Save_figures, 'true')

                 save_directory = fullfile(path_model_output, 'additionnal_plots');
                if ~exist(save_directory, 'dir')
                    mkdir(save_directory);
                end
                fig_filename = fullfile(save_directory, sprintf('friction_angle_vs_Mobility_%s.png', cleaned_headers(stat_indices(index_of_mobility_in_statistics_parameters))));
                saveas(gcf, fig_filename);
                fprintf('Figure saved: %s\n', fig_filename);
            end
        end
    catch
        disp('No mobility statistics found to plot with the average friction.');
    end


    % Average parameter module
    try
        depth_colors = lines(length(depths_for_average));
        for k = 1:length(averaged_parameters)
            % Call the get_depth_averages function to get the average values
            [time_avg, param_avg] = get_depth_averages(path_model, depths_for_average, averaged_parameters{k});
            h=figure;
            hold on;
            for i = 1:length(depths_for_average)
                plot(time_avg, param_avg{i}, 'DisplayName', sprintf('Depth = %d m', depths_for_average(i)), 'Color', depth_colors(i, :));
            end
            hold off;
            legend('Location', 'best');
            xlabel('Time (My)');
            ylabel(strrep(averaged_parameters{k}, '_', ' '));
            title(sprintf('Average %s versus Time for Depths of %s m', strrep(averaged_parameters{k}, '_', ' '), num2str(depths_for_average)));

            % Set visibility based on Display_figures
            if strcmp(Display_figures, 'true')
                set(h, 'Visible', 'on');
            else
                set(h, 'Visible', 'off');
            end

            if strcmp(Save_figures, 'true')
               save_directory = fullfile(path_model_output, 'depth_averaged_plots');
                if ~exist(save_directory, 'dir')
                    mkdir(save_directory);
                end
                % Construct a cleaned header for the specific parameter
                cleaned_header_average_parameter = regexprep(averaged_parameters{k}, '[^a-zA-Z0-9_]+', '_');
                cleaned_header_average_parameter = regexprep(cleaned_header_average_parameter, '_+$', '');  % Remove trailing underscores
                fig_filename = fullfile(save_directory, sprintf('friction_angle_vs_Mobility_%s.png', cleaned_header_average_parameter));
                saveas(gcf, fig_filename);
                fprintf('Figure saved: %s\n', fig_filename);
            end
        end
    catch
        disp('No average parameters statistics found.');
    end

    try
        if strcmp(output_all_averaged_parameters_profiles,'true')
            number_of_profiles_to_display = 10;
            %we give a dummy variable like depth to get back all the variables from
            %the averaged parameters file
            [time_var, average_params] = get_depth_averages(path_model, 'all', 'depth');
            unique_time_var = unique(time_var);
            unique_depth = unique(average_params.depth);

            % Use a colormap for better distinction between different time steps
            interval_profils = round(numel(unique_time_var)/number_of_profiles_to_display);
            vec_interval_profils = 1:interval_profils:numel(unique_time_var);

            % Loop through each variable in average_params
            variables = average_params.Properties.VariableNames;
            % Start from the third column assuming time and depth are the first two columns
            for varIndex = 3:numel(variables)
                h=figure;
                             % Set visibility based on Display_figures
                if strcmp(Display_figures, 'true')
                    set(h, 'Visible', 'on');
                else
                    set(h, 'Visible', 'off');
                end
                
                cmap = parula(numel(unique_time_var));
                for i = 1:interval_profils:numel(unique_time_var)
                    index_time_for_var= find(time_var == unique_time_var(i));

                    % Plot with distinct colors for each time step
                    if strcmp(variables(varIndex),'viscosity')
                        plot(log10(average_params.(variables{varIndex})(index_time_for_var)), ...
                            average_params.depth(index_time_for_var)./1e3, ...
                            'Color', cmap(i, :), 'LineWidth', 1.5);
                    else
                        plot(average_params.(variables{varIndex})(index_time_for_var), ...
                            average_params.depth(index_time_for_var)./1e3, ...
                            'Color', cmap(i, :), 'LineWidth', 1.5);
                    end

                    hold on;
                end

                % Reverse the y-axis for better representation of depth
                set(gca, 'YDir', 'reverse');    set(gcf, 'color', 'w');


                if strcmp(variables(varIndex),'viscosity')

                    % Add labels and title
                    xlabel('Viscosity (Pa.s)');
                    ylabel('Depth (km)');
                    title('Depth-Averaged Viscosity Over Time');
                    set(gcf, 'color', 'w');
                else
                    %             % Add labels and title for each subplot
                    %     xlabel(variables{varIndex});
                    %     ylabel('Depth (km)');
                    %     title(['Depth-Averaged ' variables{varIndex} ' Over Time']);


                    % Inside the loop where you set the xlabel and title:
                    xlabel(strrep(variables{varIndex}, '_', ' '));
                    ylabel('Depth (km)');

                    title(['Depth-Averaged ' strrep(variables{varIndex}, '_', ' ') ' Over Time']);

                end

                c = colorbar('Ticks', linspace(0, 1, number_of_profiles_to_display), ...
                    'TickLabels', arrayfun(@(x) sprintf('%.2f', x), unique_time_var(vec_interval_profils)./1e6, 'UniformOutput', false));

                % Set color bar title
                c.Label.String = 'Model Time (My)';


                %     % Add legend for better interpretation
                %     legend(arrayfun(@(x) ['Time = ' num2str(x)], unique_time_var(vec_interval_profils)./1e6, 'UniformOutput', false), 'Location', 'Best');



                if strcmp(Save_figures, 'true')
                    save_directory = fullfile(path_model_output, 'depth_averaged_plots');
                    if ~exist(save_directory, 'dir')
                        mkdir(save_directory);
                    end

                    fig_filename = fullfile(save_directory, sprintf('Depth_Averaged_%s_Over_Time.png', variables{varIndex}));
                    saveas(gcf, fig_filename);
                    fprintf('Figure saved: %s\n', fig_filename);
                end
            end
        end
    catch
        disp('Cannot plot the averaged parameters profiles');
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
            h=figure;
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

            % Set visibility based on Display_figures
            if strcmp(Display_figures, 'true')
                set(h, 'Visible', 'on');
            else
                set(h, 'Visible', 'off');
            end

            if strcmp(Save_figures, 'true')
               save_directory = fullfile(path_model_output, 'additionnal_plots');
                if ~exist(save_directory, 'dir')
                    mkdir(save_directory);
                end
                fig_filename = fullfile(save_directory, sprintf('Sinking_velocity_vs_Time.png'));
                saveas(gcf, fig_filename);
                fprintf('Figure saved: %s\n', fig_filename);
            end

            % plot periodicity for each param_avg
            h=figure;
            hold on;
            cmap = jet(numel(param_avg)); % use parula colormap with the same number of colors as the number of depths
            for i = 1:numel(param_avg)
                scatter(periods, abs(power_periods)', 100, 'o', 'filled', 'MarkerEdgeColor', 'k', 'MarkerFaceColor', cmap(i,:));
                %might be power_periods(i) to check
                xlabel('Period (years)');
                ylabel('Power (x10^3)');
                title(['Periodicity vs Time (Period = ', num2str(periods(2)), ' years)']);
                xlim([0, 600]);
                %     ylim([0, 0.08]);
            end
            legend(cellfun(@num2str, num2cell(depths_for_average./1000), 'UniformOutput', false));

            %         % Set visibility based on Display_figures
            %         if strcmp(Display_figures, 'true')
            %             set(h, 'Visible', 'on');
            %         else
            %             set(h, 'Visible', 'off');
            %         end


            %% to check later
            %         % Save the figure
            %         if strcmp(Save_figures, 'true')
            %             % Construct a cleaned header for the specific parameter
            %             cleaned_header_parameter = regexprep('Periodicity_vs_Time', '[^a-zA-Z0-9_]+', '_');
            %             cleaned_header_parameter = regexprep(cleaned_header_parameter, '_+$', '');  % Remove trailing underscores
            %
            %             fig_filename = fullfile(path_model_output, sprintf('Periodicity_vs_Time_%s.png', cleaned_header_parameter));
            %             saveas(gcf, fig_filename);
            %             fprintf('Figure saved: %s\n', fig_filename);
            %         end

        else
            %do not go throught the periodicity
        end
    catch
        disp('No periodicity statistics found.');
    end

    %Topography evolution
    if strcmp(postprocess_topography, 'true')
        try
            h=figure;hold on;
            if ~strcmp(model_geometry,'box')
            [time_elevation, elevation_topography, x_axis_interp, dip_topography] = get_topography_annulus(path_model, dt_topography, resample_topography, calculate_topography_dip, topography_smoothing_interval_for_dip_calculation);
            xlabel('Annulus Degrees [deg]')
            else
             [time_elevation, elevation_topography, x_axis_interp, dip_topography] = get_topography_box(path_model, dt_topography, resample_topography, model_length, model_height, calculate_topography_dip, topography_smoothing_interval_for_dip_calculation);
            xlabel('Model length [m]')
            end
            % Find the indices corresponding to the selected time interval
            start_index = find(time_elevation >= plot_topography_start, 1);
            end_index = find(time_elevation <= plot_topography_end, 1, 'last');
            time_elevation=time_elevation(start_index:end_index);
            elevation_topography=elevation_topography((start_index:end_index),:);
            % Plot the sorted data with adjusted x axis , here a small
            % shift in the topography of 0.2 is added so 0 appears green not blue
            surf(x_axis_interp,time_elevation,elevation_topography./1e3+0.2);shading interp;c=colorbar;demcmap('inc',[5 -8],0.1);ylabel('Time[My]');zlabel('Elevation[km]');set(gcf,'color','w');
            c.Label.String= "Elevations [km]";set(gcf,'color','w');view(2);%set(gca, 'color', 'none');grid off;set(gca,'XColor', 'none','YColor','none','ZColor','none'); % FaceLighting = 'gour

            % Set visibility based on Display_figures
            if strcmp(Display_figures, 'true')
                set(h, 'Visible', 'on');
            else
                set(h, 'Visible', 'off');
            end

            % Save the figure
            if strcmp(Save_figures, 'true')
                save_directory = fullfile(path_model_output, 'Topography');
                if ~exist(save_directory, 'dir')
                    mkdir(save_directory);
                end
                fig_filename = fullfile(save_directory, 'topography_evolution.png');
                saveas(gcf, fig_filename);
                fprintf('Figure saved: %s\n', fig_filename);
            end

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
            dip_topography=dip_topography((start_index:end_index),:);
            % Plot the sorted data with adjusted x axis
            %     For now let's use the time from statistic but this will have to be change for a resampling time, time_elevation is obsolete.
            h=figure;
            surf(x_axis_interp,time_elevation,elevation./1e3);shading interp;c=colorbar;demcmap('inc',[5 -8],0.1);ylabel('Time[My]'),xlabel('Annulus Degrees [deg]');zlabel('Elevation[km]');set(gcf,'color','w');
            c.Label.String= "Elevations [km]";set(gcf,'color','w');view(2);%set(gca, 'color', 'none');grid off;set(gca,'XColor', 'none','YColor','none','ZColor','none'); % FaceLighting = 'gour

            % Set visibility based on Display_figures
            if strcmp(Display_figures, 'true')
                set(h, 'Visible', 'on');
            else
                set(h, 'Visible', 'off');
            end

            % Save the figure
            if strcmp(Save_figures, 'true')
                                save_directory = fullfile(path_model_output, 'Dynamic_topography');
                if ~exist(save_directory, 'dir')
                    mkdir(save_directory);
                end
                fig_filename = fullfile(save_directory, 'dynamic_topography_annulus.png');
                saveas(gcf, fig_filename);
                fprintf('Figure saved: %s\n', fig_filename);
            end
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
            dip_layer=dip_layer((start_index:end_index),:);
            % Plot the sorted data with adjusted x axis
            %     For now let's use the time from statistic but this will have to be change for a resampling time, time_elevation is obsolete.
            h=figure;
            surf(x_axis_interp_layer,time_elevation_layer,elevation_layer./1e3);shading interp;c=colorbar;ylabel('Time[My]'),xlabel('Model Length [km]');zlabel('Elevation[km]');set(gcf,'color','w');
            c.Label.String= "Elevations Layer [km]";crameri('nuuk');set(gcf,'color','w');view(2);%set(gca, 'color', 'none');grid off;set(gca,'XColor', 'none','YColor','none','ZColor','none'); % FaceLighting = 'gour

            % Set visibility based on Display_figures
            if strcmp(Display_figures, 'true')
                set(h, 'Visible', 'on');
            else
                set(h, 'Visible', 'off');
            end

            % Save the figure
            if strcmp(Save_figures, 'true')

                save_directory = fullfile(path_model_output, 'Layer_topography');
                if ~exist(save_directory, 'dir')
                    mkdir(save_directory);
                end

                fig_filename = fullfile(save_directory, 'Topo_layer_evolution.png');
                saveas(gcf, fig_filename);
                fprintf('Figure saved: %s\n', fig_filename);
            end
        catch
            disp('Topography layer files not found.');
        end
    end

    %Calculated dip of topography or layer evolution
    if strcmp(calculate_topography_dip, 'true')
        try
            %       Plot the dip over time
            x_axis_dip = x_axis_interp(2:end);
            h=figure;hold on;
            if ~strcmp(model_geometry,'box')
            xlabel('Annulus Degrees [deg]')
            else
            xlabel('Model length [m]')
            end
            surf(x_axis_dip, time_elevation, dip_topography);
            shading interp;
            c = colorbar;
            % Find the max value of the dip
            max_value = max(max(abs(dip_topography)));
            % Set color limits based on the maximum dip
            clim([-max_value, max_value]);
            crameri('vik');
            ylabel('Time[My]');
            zlabel('Dip Topography (°)');
            set(gcf, 'color', 'w');
            c.Label.String = "Dip Topography (°)";
            view(2);

            disp('Dip topography calculation ');

            % Set visibility based on Display_figures
            if strcmp(Display_figures, 'true')
                set(h, 'Visible', 'on');
            else
                set(h, 'Visible', 'off');
            end

            % Save the figure
            if strcmp(Save_figures, 'true')
                          save_directory = fullfile(path_model_output, 'Topography');
                if ~exist(save_directory, 'dir')
                    mkdir(save_directory);
                end
                fig_filename = fullfile(save_directory, 'Topography_dip_evolution.png');
                saveas(gcf, fig_filename);
                fprintf('Figure saved: %s\n', fig_filename);
            end

        catch
            disp('Dip topography calculation not asked or not possible.');
        end
    elseif strcmp(calculate_topography_layer_dip, 'true')
        try
            %       Plot the dip over time
            x_axis_dip_layer = x_axis_interp_layer(2:end);
            h=figure;
            surf(x_axis_dip_layer, time_elevation_layer, dip_layer);
            shading interp;
            c = colorbar;
            % Find the max value of the dip
            max_value = max(max(abs(dip_layer)));
            % Set color limits based on the maximum dip
            clim([-max_value, max_value]);
            crameri('vik');
            ylabel('Time[My]'), xlabel('Annulus Degrees [deg]');
            zlabel('Dip Layer');
            set(gcf, 'color', 'w');
            c.Label.String = "Dip Layer";
            view(2);
            disp('Dip layer calculation ');

            % Set visibility based on Display_figures
            if strcmp(Display_figures, 'true')
                set(h, 'Visible', 'on');
            else
                set(h, 'Visible', 'off');
            end

            % Save the figure
            if strcmp(Save_figures, 'true')
                save_directory = fullfile(path_model_output, 'Layer_topography');
                if ~exist(save_directory, 'dir')
                    mkdir(save_directory);
                end
                fig_filename = fullfile(save_directory, 'Topography_layer_evolution.png');
                saveas(gcf, fig_filename);
                fprintf('Figure saved: %s\n', fig_filename);
            end

        catch
            disp('Dip layer calculation not not asked or not possible.');
        end
        try
            %       Plot the dip evolution over time as a mean and standart deviation
            %       Calculate the mean and standard deviation of dip values over time

            %        Only produce statistic of slope when significant change of slope occurs
            %        the user can set a slope_gradient_threshold (e.g 0.1 equivalent to 10prct),
            %        the statistic will be produced for slope having a gradient higher than this number
            % Define the slope gradient threshold
            slope_gradient_threshold = 0.02;
            minimum_slope_threshold = 10;

            %       filter with slope
            slope_filter = abs(dip_layer)>minimum_slope_threshold;
            % Calculate the absolute gradient of dip_layer
            gradient_filter = abs(gradient(dip_layer))>slope_gradient_threshold;


            % Create a figure showing the slope areas used for the plot
            h=figure;
            surf(x_axis_dip_layer, time_elevation_layer,slope_filter.*gradient_filter);
            shading flat;
            colorbar;
            clim([0, 1]);
            view(2);
            % Add labels and a title
            xlabel('Model length [km]');
            ylabel('Time [yr]');
            title('Slope Areas Used for Plot');

            % Set the figure background color to white
            set(gcf, 'color', 'w');

            % Optionally, you can add grid lines
            grid on;

            % Set visibility based on Display_figures
            if strcmp(Display_figures, 'true')
                set(h, 'Visible', 'on');
            else
                set(h, 'Visible', 'off');
            end

            % Save the figure
            if strcmp(Save_figures, 'true')
                save_directory = fullfile(path_model_output, 'Topography');
                if ~exist(save_directory, 'dir')
                    mkdir(save_directory);
                end
                fig_filename = fullfile(save_directory, 'topo_layer_threshold_areas.png');
                saveas(gcf, fig_filename);
                fprintf('Figure saved: %s\n', fig_filename);
            end

            dip_layer_tracked = slope_filter.*gradient_filter;
            dip_layer_with_gradient_threshold = dip_layer.*dip_layer_tracked;
            dip_layer_with_gradient_threshold(dip_layer_with_gradient_threshold == 0) = NaN;

            %         Create an histogram of the dip values taken into account the filters
            %         h=figure;
            %         % Specify the number of bins for the histogram
            %         num_bins = 20;
            %
            %         % Create a histogram
            %         histogram(abs(dip_layer_with_gradient_threshold), num_bins);
            %
            %         % Add labels and a title
            %         xlabel('Dip Layer Values');
            %         ylabel('Frequency');
            %         title('Histogram of Dip Layer Values');


            mean_dip = nanmean(abs(dip_layer_with_gradient_threshold), 2);
            std_dip = nanstd(abs(dip_layer_with_gradient_threshold), 0, 2);
            non_nan_values = sum(~isnan(dip_layer_with_gradient_threshold), 2);


            % Your mean and std_dip vectors
            y = mean_dip'; % Replace with your mean dip values
            std_dev = std_dip'; % Replace with your std dip values


            %         % Calculate upper and lower curves for shading without confidency interval
            %         curve1 = y + std_dev;
            %         curve2 = y - std_dev;

            % Calculate upper and lower curves with confidency interval

            %         In case of no dip detetected replace nan values by 0
            non_nan_values(isnan(non_nan_values))=0;
            for i = 1:numel(non_nan_values)
                SEM = std_dev(i) / sqrt(non_nan_values(i)); % Standard Error

                % Calculate ts (T-Score) for the current non_nan_values
                %             confidence interval is set at 70%
                ts = tinv([0.15 0.85], non_nan_values(i) - 1);

                % Calculate Confidence Intervals using ts
                curve1(i) = y(i) + ts(2) * SEM; % Upper bound of Confidence Interval
                curve2(i) = y(i) + ts(1) * SEM; % Lower bound of Confidence Interval
            end


            % Create x2 and inBetween vectors for shading
            x2 = [time_elevation_layer, fliplr(time_elevation_layer)];
            inBetween = [curve1, fliplr(curve2)];

            h=figure;
            % Plot the shaded area
            fill(x2, inBetween, 'b', 'FaceAlpha', 0.3); % Adjust color and transparency as needed

            hold on;

            % Plot the mean dip values
            plot(time_elevation_layer, mean_dip, 'r-', 'LineWidth', 2); % Adjust color and line style as needed

            xlabel('Time');
            ylabel('Mean Dip °');
            title('Mean Dip with Standard Deviation in Shaded Area');
            legend('Standard Deviation', 'Mean Dip', 'Location', 'Best');
            set(gcf, 'color', 'w');
            grid on;
            disp('Dip mean and standart deviation calculation')

            % Set visibility based on Display_figures
            if strcmp(Display_figures, 'true')
                set(h, 'Visible', 'on');
            else
                set(h, 'Visible', 'off');
            end

            % Save the figure
            if strcmp(Save_figures, 'true')
                save_directory = fullfile(path_model_output, 'Layer_topography');
                if ~exist(save_directory, 'dir')
                    mkdir(save_directory);
                end
                fig_filename = fullfile(save_directory, 'topo_layer_mean_dip_versus_time.png');
                saveas(gcf, fig_filename);
                fprintf('Figure saved: %s\n', fig_filename);
            end

        catch
            disp('Dip mean and standart deviation calculation not possible.');
        end
    end



    % Surface heatflux evolution
    if strcmp(postprocess_heatflux, 'true')
        try
            [time_heatflux, heatfluxmap, x_axis_interp_heatflux] = get_heatflux_annulus(path_model, dt_heatflux, resample_heatflux);
            % Plot the sorted data with adjusted x axis
            h=figure;
            surf(x_axis_interp_heatflux,time(1:time_heatflux:size(elevation,1)),heatfluxmap.*1e3);shading interp;c = colorbar;crameri('lajolla',12);ylabel('Time[My]'),xlabel('Annulus Degrees [deg]');zlabel('Surface heat flux[W]');caxis([0 120]);set(gcf,'color','w');
            c.Label.String= "Surface heat flux [mW/m]";set(gcf,'color','w');view(2);%set(gca, 'color', 'none');grid off;set(gca,'XColor', 'none','YColor','none','ZColor','none'); % FaceLighting = 'gour
            % Save the figure

            % Set visibility based on Display_figures
            if strcmp(Display_figures, 'true')
                set(h, 'Visible', 'on');
            else
                set(h, 'Visible', 'off');
            end

            if strcmp(Save_figures, 'true')
                save_directory = fullfile(path_model_output, 'Heat_flux');
                if ~exist(save_directory, 'dir')
                    mkdir(save_directory);
                end
                fig_filename = fullfile(save_directory, 'heat_flux_annulus.png');
                saveas(gcf, fig_filename);
                fprintf('Figure saved: %s\n', fig_filename);
            end

        catch
            disp('Heatflux files not found. Skipping heatflux.');
        end
    end

    %Plot subduction and Plumes statistics
    if strcmp(run_spherical_additional_postprocess, 'true')
        filename_geofeatures = ['Geofeatures_' model_title '.csv'];
        fullFilePath_geofeatures = fullfile(path_model_output, filename_geofeatures);
        % Check if the file exists
        file_geofeatures_exists = exist(fullFilePath_geofeatures, 'file');
        init_step=0;
        if strcmp(restart_additional_postprocess,'auto')
            [init_step, max_step] = get_last_timestep(path_model, path_model_output,output_additional_maps_figures);
            if init_step == max_step
                disp('The model has already finished computing the spherical additional postprocess.');
            else
                disp(['Additional spherical postprocessing will start or restart from init_step ' num2str(init_step)]);

                % If it doesn't exist then run the additional postprocess function
                if  strcmp(restart_additional_postprocess,'auto') && (init_step < max_step)
                    [geofeatures] = get_geodynamics_features_statistics(path_model_input, additional_postprocesses,plumes_depths_tracking,subduction_depths_tracking,...
                        plumes_non_adiabatic_tracking_temperature,subduction_non_adiabatic_tracking_temperature,trenches_elevation_threshold,path_model_output,...
                        compositional_field_name_of_continents,output_figures_for_spherical_additional_postprocess,output_additional_maps_figures,write_geofeatures_statistics,...
                        remove_subductions_to_oceanic_age,init_step,Interval_of_time_output_for_additional_postprocess,topography_correction,reference_time_My_to_Ma,plot_continents_border_from_reconstruction,...
                        additional_fields_to_load,additional_fields_threshold,default_threshold_fields,additional_fields_depths_to_visualize);
                    % If the file exists then plot 'subduction_and_plume_statistics'
                end
            end
        end

        try
            if file_geofeatures_exists ~= 0
                fprintf('The file %s was found. It will plotted.\n', fullFilePath_geofeatures);
                %Lets remove subduction and plumes statistics from the additional
                %postprocess if the file already exist
                %obsolete  : additional_postprocesses(strcmp(additional_postprocesses, 'subduction_and_plume_statistics')) = [];
                geofeatures = import_geofeature_statistics(fullFilePath_geofeatures);
            end
            % Get unique depths
            uniqueDepths = unique(geofeatures.Depths);
            if all(isnan(uniqueDepths))
                iteration_over_deth = 1;
            else
                iteration_over_depth = length(uniqueDepths);
            end
            % Create a cell array to store tables
            depthTables = cell(size(uniqueDepths));
            % Loop over each depth and create a table
            for i = 1:iteration_over_depth
                %%
                depth = uniqueDepths(i);
                depthTables{i} = geofeatures(geofeatures.Depths == depth, :);
                depthTables{i}.Depths;
                % Plot the data for the first depth
                h=figure;

                % Create the first y-axis for 'Number of Subductions' and 'Number of Subduction Connected'
                yyaxis left;ylim([0, 150]);
                % Plot the first line with solid line ('k-')
                plot(depthTables{i}.Time, depthTables{i}.NumSubductions, 'k-', 'LineWidth', 2);
                hold on;  % Hold the plot to add another line
                % Plot the second line with dashed line ('k--')
                plot(depthTables{i}.Time, depthTables{i}.NumSubduction_connected, 'k-*', 'LineWidth', 2);
                hold on;
                plot(depthTables{i}.Time, depthTables{i}.Maxferet_mean./1e2,'b-', 'LineWidth', 2);
                % yl1 = ax1.YLim;
                % ax1.YTick = linspace(yl1(1), yl1(2), length(ax1.YTick));
                hold on;

                % Create the second y-axis for 'Number of Plumes'
                yyaxis right;ylim([0, 250]);
                % Plot the third line with a different color and line style ('r-')
                plot(depthTables{i}.Time, depthTables{i}.NumPlumes, 'r-', 'LineWidth', 2);
                % Add labels and title
                xlabel('Time (Myr)');
                yyaxis left; ylabel('Number of Subductions and Mean maximum feret diamater (.1e2 [km])');
                yyaxis right; ylabel('Number of Plumes');
                title(['Number of Subductions and Plumes Over Time for Depth ' num2str(uniqueDepths(i)) ' km']);

                % Add grid and legend
                grid on;
                set(gcf, 'color', 'w');
                legend({'Number of Subductions', 'Number of Subduction Connected','Mean maximum feret diamater', 'Number of Plumes'});

                % Adjust the appearance of the plot
                set(gca, 'FontSize', 12);
                %%

                % Set visibility based on Display_figures
                if strcmp(Display_figures, 'true')
                    set(h, 'Visible', 'on');
                else
                    set(h, 'Visible', 'off');
                end

                % Save the figure
                if strcmp(Save_figures, 'true')
                   save_directory = fullfile(path_model_output, 'additionnal_paraview_plots');
                    if ~exist(save_directory, 'dir')
                        mkdir(save_directory);
                    end
                    fig_filename = fullfile(save_directory, sprintf('Subduction_plumes_%d.png', depth));
                    saveas(gcf, fig_filename);
                    fprintf('Figure saved: %s\n', fig_filename);
                end

            end
        catch
            disp('Subduction and plumes numbers could not or was not plotted.');
        end
        if any(strcmp(additional_postprocesses, 'continents_VRMS'))
            try
                %Plot the Vrms for the surface and continents only
                % Get unique depths
                [uniqueTime,index_uniqueTime] = unique(geofeatures.Time);
                Vrms_surface = geofeatures.V_rms_surface(index_uniqueTime);
                Vrms_continents = geofeatures.V_rms_continents(index_uniqueTime);

                % Plot the data for the first depth
                h=figure;
                plot(uniqueTime, Vrms_continents, 'b-', 'LineWidth', 2, 'DisplayName', 'VRMS Continents');
                hold on;
                plot(uniqueTime, Vrms_surface, 'r-', 'LineWidth', 2, 'DisplayName', 'VRMS Surface');

                % Adding labels and title
                xlabel('Time (My)');
                ylabel('VRMS (m.yr^-^1)');
                title('VRMS over Time for Surface and Continents only');
                legend('show');
                grid on;
                set(gcf, 'color', 'w');

                % Set visibility based on Display_figures
                if strcmp(Display_figures, 'true')
                    set(h, 'Visible', 'on');
                else
                    set(h, 'Visible', 'off');
                end

                if strcmp(Save_figures, 'true')
                     save_directory = fullfile(path_model_output, 'additionnal_paraview_plots');
                    if ~exist(save_directory, 'dir')
                        mkdir(save_directory);
                    end
                    fig_filename = fullfile(save_directory, 'Vrms_continents_and_surface.png');
                    saveas(gcf, fig_filename);
                    fprintf('Figure saved: %s\n', fig_filename);
                end

                try
                    %Recalculate the mobility with the Vrms of the continents only
                    index_of_RMS_velocity_in_statistics_parameters = find(contains(statistic_parameters, 'RMS velocity'),1);
                    if isempty(index_of_RMS_velocity_in_statistics_parameters)
                        disp('Could not calculate the mobility of continents, RMS velocity has to be added to statistic_parameters')
                    end

                    % Vrms global at Vrms continents time resolution
                    %we need to remove the initial refinement steps from time from statistics
                    [uniquetime_stats,index_uniquetime_stats] = unique(time);
                    Vrms_global = data_stats(:,stat_indices(index_of_RMS_velocity_in_statistics_parameters));
                    Vrms_less_IRS = Vrms_global(index_uniquetime_stats);
                    V_RMS_resol_mobility= interp1(time(index_uniquetime_stats)./1e6,Vrms_less_IRS,uniqueTime);
                    Continents_mobility =Vrms_continents./V_RMS_resol_mobility;

                    index_of_mobility_in_statistics_parameters = find(contains(statistic_parameters, 'Mobility'),1);
                    if isempty(index_of_mobility_in_statistics_parameters)
                        disp('Could not calculate the mobility of continents, "Mobility" has to be added to statistic_parameters')
                    end

                    Mobility_global = data_stats(:,stat_indices(index_of_mobility_in_statistics_parameters));

                    Mobility_less_IRS = Mobility_global(index_uniquetime_stats);

                    Mobility_resol_mobility= interp1(time(index_uniquetime_stats)./1e6,Mobility_less_IRS,uniqueTime);

                    %Plot the continents mobility
                    h=figure;
                    plot(uniqueTime, Continents_mobility, 'b-', 'LineWidth', 2, 'DisplayName', 'Continents Mobility');
                    hold on;
                    plot(uniqueTime,Mobility_resol_mobility, 'r-', 'LineWidth', 2, 'DisplayName', 'Full Mobility');
                    xlim([min(uniqueTime), max(uniqueTime)]);
                    % Adding labels and title
                    xlabel('Time (My)');
                    ylabel('Continents Mobility');
                    title('Continents Mobility over Time');
                    legend('show');
                    grid on;
                    set(gcf, 'color', 'w');

                    % Set visibility based on Display_figures
                    if strcmp(Display_figures, 'true')
                        set(h, 'Visible', 'on');
                    else
                        set(h, 'Visible', 'off');
                    end

                    if strcmp(Save_figures, 'true')
                     save_directory = fullfile(path_model_output, 'additionnal_paraview_plots');
                    if ~exist(save_directory, 'dir')
                        mkdir(save_directory);
                    end
                        fig_filename = fullfile(save_directory, 'Continents_Mobility.png');
                        saveas(gcf, fig_filename);
                        fprintf('Figure saved: %s\n', fig_filename);
                    end

                catch
                    disp('Could not calculate and plot the continents mobility.');

                end
            catch
                disp('Could not plot the Vrms of the continents.');

            end

        if any(strcmp(additional_postprocesses, 'oceanic_age_statistics'))
            try
               filename_oceans = ['Geo_oceans_' model_title '.csv'];
               fullFilePath_geooceans = fullfile(path_model_output, filename_oceans);
               oceanic_age_density = import_layers_sphere(fullFilePath_geooceans);
               oceanic_age_density_reshape = reshape(oceanic_age_density.Time_with_Pdf_oceanic_age,[401,numel(oceanic_age_density.Time_with_Pdf_oceanic_age)/401]);
               oceanic_age_density_time = oceanic_age_density_reshape(1,:);
               oceanic_age_density_reshape = oceanic_age_density_reshape(2:end, :);
               oceanic_age_distribution = 1:400;
                h = figure;
                pcolor(oceanic_age_density_time, oceanic_age_distribution, oceanic_age_density_reshape);
                shading flat;
                
                % Add labels and title
                xlabel('Time (Myr)');
                ylabel('Oceanic Age Distribution (Ma)');
                title('Oceanic Age Density Distribution with Kernel:  epanechnikov, bw = 0.3 ');
                
                % Add colorbar with better readability
                c = colorbar;
                c.Label.String = 'Density distribution';
                cmap = crameri('-roma');
                colormap(cmap);
                caxis([0, 0.01]); % Set colorbar limits
                set(gcf, 'color', 'w');
               
                % Set visibility based on Display_figures
                if strcmp(Display_figures, 'true')
                    set(h, 'Visible', 'on');
                else
                    set(h, 'Visible', 'off');
                end
        
                if strcmp(Save_figures, 'true')
                    save_directory = fullfile(path_model_output, 'additionnal_paraview_plots');
                    if ~exist(save_directory, 'dir')
                        mkdir(save_directory);
                    end
                    fig_filename = fullfile(save_directory, 'Continents_Mobility.png');
                    saveas(gcf, fig_filename);
                    fprintf('Figure saved: %s\n', fig_filename);
                end

            catch
                disp('Could no plot the Oceanic age density distribution.');
            end
        end


        end
    end
end


if strcmp(Display_figures, 'true')
    % Tile all open figures
    tilefigs;
end
