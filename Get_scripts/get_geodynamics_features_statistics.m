function [geofeatures] = get_geodynamics_features_statistics(path_model_input, additional_postprocesses,plumes_depths_tracking,subduction_depths_tracking,plumes_non_adiabatic_tracking_temperature,...
    subduction_non_adiabatic_tracking_temperature,trenches_elevation_threshold,path_model_output,compositional_field_name_of_continents,output_figures_for_spherical_additional_postprocess,output_additional_maps_figures,...
    write_geofeatures_statistics,remove_subductions_to_oceanic_age,init_step,Interval_of_time_output_for_additional_postprocess,topography_correction,reference_time_My_to_Ma,plot_continents_border_from_reconstruction,...
    additional_fields_to_load,additional_fields_threshold,default_threshold_fields,additional_fields_depths_to_visualize,visualize_model_at_specific_time,end_step,Display_figures,plot_an_additional_field_on_top_of_geofeatures)

warning('off','all');

allowed_parameters = {'subduction_and_plume_statistics','oceanic_age_statistics','continents_VRMS','melt_statistics'};
% Check if averaged_parameter is allowed
if ~all(ismember(additional_postprocesses, allowed_parameters))
    error('Invalid parameter selected. The allowed parameters are: %s', strjoin(allowed_parameters, ', '));
end


allowed_maps = {'topography','geofeatures','strain_rate','oceanic_age','divergence'};
% Check if averaged_parameter is allowed
if ~all(ismember(output_additional_maps_figures, allowed_maps))
    error('Invalid maps selected. The allowed maps are: %s', strjoin(allowed_maps, ', '));
end

geofeatures = table();

% Split the path into parts
parts = strsplit(path_model_input, '/');

% Extract the last part
model_title = parts{end-1};

%initiate these variable for further checking
num_files_surface=[];
num_files_lithosphere=[];
num_files_depths=[];
% Default values for uninitialized variables
default_value = NaN;

% Initialize the variables
numPlumes = default_value;
numSubductions = default_value;
numSubduction_connected = default_value;
current_depth = default_value;
Maxferet_mean = default_value;
Maxferet_max = default_value;
Maxferet_min = default_value;
nferet = default_value;
V_rms_continents=default_value;
V_rms_surface=default_value;
files_required_1 =[];
files_required_2 =[];
files_required_3 =[];

RMS_melt_surface = default_value;
RMS_melt_depths = default_value;
RMS_melt_sublithosphere = default_value;

surface_additional_fields_to_load_map = [];
lithosphere_additional_fields_to_load_map = [];
depths_additional_fields_to_load_map = [];
surface_additional_fields = [];
lithosphere_additional_fields = [];
depths_additional_fields = [];

additional_fields_depths_to_visualize_ite = additional_fields_depths_to_visualize;
additional_fields_threshold_ite = additional_fields_threshold;


%Variable for velocity vectors
scaling_velocity_factor_surface=0.2;
scaling_velocity_factor_depths=0.01;
%Variable colors for plot
CMap_plume = [1 0 0 ;  1 1 1];
CMap_subduction = [0 0 1 ;  1 1 1];
CMap_boundaries = [ 0.0471 0.4392 0.2275;  1 1 1];
CMap_continents = [ 0.4745 0.3647 0.3020;  1 1 1];
CMap_continents_dark = [ 0 0 0;  1 1 1];
CMap_additional_fields = [1 0.2 0 ;  1 1 1];
CMap_additional_fields_for_geofeatures = [0.8 0 0.8 ;  1 1 1];


%Required files for each scheme
if any(strcmp(additional_postprocesses, 'subduction_and_plume_statistics'))
    files_required_1 = {'surface','depths'};
end
if any(strcmp(additional_postprocesses, 'oceanic_age_statistics')) || any(strcmp(output_additional_maps_figures, 'oceanic_age')) || any(strcmp(output_additional_maps_figures, 'divergence'))
    files_required_2 = {'lithosphere'};
end
if any(strcmp(additional_postprocesses, 'continents_VRMS')) ||any(strcmp(output_additional_maps_figures, 'topography')) || strcmp(remove_subductions_to_oceanic_age, 'true')
    files_required_3 = {'surface'};
end
total_files_required = [files_required_1 files_required_2 files_required_3];
files_requirement = unique(total_files_required);

% Check if subduction_and_plume_number is true
if any(strcmp(files_requirement, 'surface'))
    file_surface = dir(fullfile(path_model_input, 'surface_*'));
    filenames_surface = {file_surface.name};
    num_files_surface = numel(filenames_surface);
    disp(['Found ' num2str(num_files_surface) ' surface files.']);

end

if any(strcmp(files_requirement, 'depths'))
    file_depths = dir(fullfile(path_model_input, 'depths_*'));
    filenames_depths = {file_depths.name};
    num_files_depths = numel(filenames_depths);
    disp(['Found ' num2str(num_files_depths) ' depths files.']);

end


if any(strcmp(files_requirement, 'lithosphere'))
    file_lithosphere = dir(fullfile(path_model_input, 'lithosphere_*'));
    filenames_lithosphere = {file_lithosphere.name};
    num_files_lithosphere = numel(filenames_lithosphere);
    % If no files are there write a message, alternatively, wrrite how many
    % files you found

    % Process lithosphere files
    disp(['Found ' num2str(num_files_lithosphere) ' lithosphere files for the additional postprocess.']);
    % Add code to handle the presence of lithosphere files if needed
end

% To optimize the code, later, we could use some sorting
% file_counts = [num_files_surface, num_files_lithosphere, num_files_depths];
%
% % Sort indices based on the number of files in descending order
% [sorted_file_counts, sorted_file_indices] = sort(file_counts, 'descend');

if any(strcmp(files_requirement, 'surface'))
    if (num_files_surface == 0)
        disp('No relevant surface files found. Exiting the function.');
        return;
    end
end


if any(strcmp(files_requirement, 'lithosphere'))
    if (num_files_lithosphere == 0)
        disp('No relevant lithosphere files found. Exiting the function.');
        return;
    end
end

if any(strcmp(files_requirement, 'depths'))
    if (num_files_depths == 0)
        disp('No relevant depths surface found. Exiting the function.');
        return;
    end
end

%check how many iteration we should do, but do not overshoot if we dont
%have the same number of file between surface lithosphere and depths for
%some of the plots
min_num_files = min([num_files_lithosphere, num_files_surface, num_files_depths]);
max_num_files = max([num_files_lithosphere, num_files_surface, num_files_depths]);


% Check if the number of files for surface and depths is different
if num_files_surface ~= num_files_depths
    fprintf('Number of surface files (%d) is different from the number of depth files (%d), the satistics will only be counted until file number (%d).', num_files_surface, num_files_depths,min_num_files);
end

count_ite = 0;
count_ite_with_steps=0;
Reference_mean_correction_steps=2;

% Create a grid of longitude and latitude values
[Xeq, Yeq] = meshgrid(-180:1:180, -90:1:90);
% and a lower resolutin for velocity glyphs
[Xeq_lr, Yeq_lr] = meshgrid(-180:5:180, -90:5:90);

if any(strcmp(files_requirement, 'surface'))

    %We want the initial position of the continent to be plotted on top of the maps
    surface_init = import_layers_sphere(fullfile(path_model_input, filenames_surface{Reference_mean_correction_steps}));
    % Sort the data in the desired order
    x_surface_init = surface_init.("Points:0");
    y_surface_init = surface_init.("Points:1");
    z_surface_init = surface_init.("Points:2");

    % Calculate the radial distance from the center of the Earth
    r_surface_init = sqrt(x_surface_init.^2 + y_surface_init.^2 + z_surface_init.^2);

    surface_topography_init = r_surface_init-6371e3;

    % Calculate the azimuthal angle
    theta_surface_init = atan2(y_surface_init, x_surface_init);
    [theta_sorted_init, theta_index_init] = sort(theta_surface_init);
    x_sorted_surface_init = x_surface_init(theta_index_init);
    y_sorted_surface_init = y_surface_init(theta_index_init);
    z_sorted_surface_init = z_surface_init(theta_index_init);
    surface_init_topography_sorted = surface_topography_init(theta_index_init);
    %     surface_init_continent_sorted = surface_init.continent(theta_index_init);

    % Check if the field exists and assign the parameter accordingly
    if ismember(compositional_field_name_of_continents, surface_init.Properties.VariableNames)
        surface_init_continent_sorted = surface_init.(compositional_field_name_of_continents)(theta_index_init);
    else
        % Handle the case where the field does not exist
        surface_init_continent_sorted = 'NaN';  % or assign a default value
        % Display a message
        disp('Warning: There is no field assigned to continents. This is fine if your model has no continents, otherwise please assign a field so they can be tracked and projected.');
    end

    % Convert Cartesian coordinates to longitude and latitude
    longitude_surface_init = atan2(y_sorted_surface_init, x_sorted_surface_init) * 180 / pi;
    latitude_surface_init = asin(z_sorted_surface_init ./ sqrt(x_sorted_surface_init.^2 + y_sorted_surface_init.^2 + z_sorted_surface_init.^2)) * 180 / pi;

    if ~strcmp(surface_init_continent_sorted,'NaN')
        contour_continents_init = griddata(longitude_surface_init, latitude_surface_init, surface_init_continent_sorted, Xeq, Yeq);
    else
        %assign continent is null everywhere
        contour_continents_init = zeros(size(Xeq));
    end

    %Will be use to correct the topography (I should make it an option)
    mean_surface_init_topography = mean(surface_init_topography_sorted);
end

% For restart
if isempty(visualize_model_at_specific_time)
    if init_step~=0
        start_step = init_step;
         final_step = max_num_files;
    else
        start_step = 1;
         final_step = max_num_files;
    end
elseif ~isempty(visualize_model_at_specific_time)
    start_step = end_step;
    final_step = end_step;

else 
    start_step = 1;
    final_step = max_num_files;
end

for i = start_step:final_step

    count_ite = count_ite+1;
    count_ite_with_steps = start_step+count_ite+1;
    fprintf('Processing additional statistics files %d of %d\n', i, min_num_files)
    fprintf('Processing additional maps files %d of %d\n', i, max_num_files)

    file_index = i;
    try
        % Import only the necessary data
        if  any(strcmp(files_requirement, 'surface')) && i<=num_files_surface
            surface = import_layers_sphere(fullfile(path_model_input, filenames_surface{file_index}));
        end
        if  any(strcmp(files_requirement, 'depths')) && i<=num_files_depths
            depths = import_layers_sphere(fullfile(path_model_input, filenames_depths{file_index}));
        end
        if any(strcmp(files_requirement, 'lithosphere'))&& i<=num_files_lithosphere
            lithosphere = import_layers_sphere(fullfile(path_model_input, filenames_lithosphere{file_index}));
        end

        % Lets define where the time should be read
        if any(strcmp(files_requirement,'surface'))
            if isfield(surface, 'Time')
                time_str = sprintf('Time: %.2f Myr', surface.Time(1) / 1e6);
                time_for_output_structures = surface.Time(1) / 1e6;
            else
                % If the field time doesn't exist
                time_for_output_structures = (i-1) * Interval_of_time_output_for_additional_postprocess;
                time_str = sprintf('Time: %.2f Myr', time_for_output_structures);
            end
        elseif any(strcmp(files_requirement,'lithosphere'))
            if isfield(lithosphere, 'Time')
                time_str = sprintf('Time: %.2f Myr', lithosphere.Time(1) / 1e6);
                time_for_output_structures = lithosphere.Time(1) / 1e6;
            else
                % Set a default time value for lithosphere
                time_for_output_structures = (i-1) * Interval_of_time_output_for_additional_postprocess; % Adjust this default value if needed
                time_str = sprintf('Time: %.2f Myr', time_for_output_structures);
            end
        elseif any(strcmp(files_requirement,'depths'))
            if isfield(depths, 'Time')
                time_str = sprintf('Time: %.2f Myr', depths.Time(1) / 1e6);
                time_for_output_structures = depths.Time(1) / 1e6;
            else
                % Set a default time value for depths
                time_for_output_structures = (i-1) * Interval_of_time_output_for_additional_postprocess; % Adjust this default value if needed
                time_str = sprintf('Time: %.2f Myr', time_for_output_structures);
            end
        end


        if strcmp(plot_continents_border_from_reconstruction,'true')
            % Convert model time to corresponding filename time
            filename_time = reference_time_My_to_Ma - time_for_output_structures;

            % Round down to the nearest integer to find the closest matching filename
            rounded_filename_time = floor(filename_time);

            % Create the filename string
            filename_time_str = sprintf('../data/reconstructed_%.2fMa.shp', rounded_filename_time);

            % Load the shapefile if needed
            shape_continents = shaperead(filename_time_str,'UseGeoCoords',true);
        end

        % Define the downsampling factor
        downsample_factor = 1;

        if  any(strcmp(files_requirement, 'surface'))

            % Sort the data in the desired order
            x_surface = surface.("Points:0");
            y_surface = surface.("Points:1");
            z_surface = surface.("Points:2");

            % Calculate the radial distance from the center of the Earth
            r_surface = sqrt(x_surface.^2 + y_surface.^2 + z_surface.^2);

            surface_topography = r_surface-6371e3;

            % Calculate the azimuthal angle
            theta_surface = atan2(y_surface, x_surface);

            % Sort the data based on the azimuthal angle
            [theta_sorted, theta_index] = sort(theta_surface);
            x_sorted_surface = x_surface(theta_index);
            y_sorted_surface = y_surface(theta_index);
            z_sorted_surface = z_surface(theta_index);
            r_surface_sorted = r_surface(theta_index);
            surface_velocity0_sorted = surface.("velocity:0")(theta_index);
            surface_velocity1_sorted = surface.("velocity:1")(theta_index);
            surface_velocity2_sorted = surface.("velocity:2")(theta_index);
            surface_topography_sorted = surface_topography(theta_index);
            surface_strain_rate_sorted = surface.strain_rate(theta_index);

            if ~strcmp(additional_fields_to_load{1},'')
                additional_fields_to_load_map = cell(numel(additional_fields_to_load));
                for f = 1:numel(additional_fields_to_load)
                    if ismember(additional_fields_to_load{f}, surface.Properties.VariableNames)
                        % Retrieve the data using a temporary variable
                        temp_data = surface.(additional_fields_to_load{f})(theta_index);
                        % Assign to the cell array using curly braces
                        surface_additional_fields_to_load_sorted{f} = temp_data;
                    end
                end
            end


            % Check if the field exists and assign the parameter accordingly
            if ismember(compositional_field_name_of_continents, surface.Properties.VariableNames)
                surface_continent_sorted = surface.(compositional_field_name_of_continents)(theta_index);
            else
                % Handle the case where the field does not exist
                surface_continent_sorted = 'NaN';  % or assign a default value
            end

            % Check if 'continents_VRMS' was selected and the field is NaN
            if ismember('continents_VRMS', additional_postprocesses) && strcmp(surface_continent_sorted,'NaN')
                % Display a message
                disp('Warning: There is currently no field assigned as continents to calculate the VRMS of continents!');
            end

            final_surface_topography_sorted=surface_topography_sorted;


            x_surface_downsampled = x_sorted_surface(1:downsample_factor:end);
            y_surface_downsampled = y_sorted_surface(1:downsample_factor:end);
            z_surface_downsampled = z_sorted_surface(1:downsample_factor:end);
            surface_topography_surface_downsampled = final_surface_topography_sorted(1:downsample_factor:end);
            % Convert Cartesian coordinates to longitude and latitude
            longitude_surface_downsampled = atan2(y_surface_downsampled, x_surface_downsampled) * 180 / pi;
            latitude_surface_downsampled = asin(z_surface_downsampled ./ sqrt(x_surface_downsampled.^2 + y_surface_downsampled.^2 + z_surface_downsampled.^2)) * 180 / pi;
        end

        if any(strcmp(files_requirement, 'depths'))  && i<=num_files_depths
            x_depths = depths.("Points:0");
            y_depths = depths.("Points:1");
            z_depths = depths.("Points:2");

            theta_depths = atan2(y_depths, x_depths);
            [theta_sorted_depths, theta_index_depths] = sort(theta_depths);
            x_sorted_depths = x_depths(theta_index_depths);
            y_sorted_depths = y_depths(theta_index_depths);
            z_sorted_depths = z_depths(theta_index_depths);
            depths_sorted_depths = depths.depth(theta_index_depths);
            unique_depths = unique(depths_sorted_depths);
            depths_nonadiabatic_temperature_sorted = depths.nonadiabatic_temperature(theta_index_depths);
            %             additional_fields_to_load = {'llsvps','llsvps'};
            if ~strcmp(additional_fields_to_load{1},'')
                additional_fields_to_load_map = cell(numel(unique_depths), numel(additional_fields_to_load));
                for f = 1:numel(additional_fields_to_load)
                    if ismember(additional_fields_to_load{f}, depths.Properties.VariableNames)
                        % Retrieve the data using a temporary variable
                        temp_data = depths.(additional_fields_to_load{f})(theta_index_depths);
                        % Assign to the cell array using curly braces
                        additional_fields_to_load_sorted{f} = temp_data;
                    end
                end
            end


            for u = 1:numel(unique_depths)
                idx_depths = find(depths_sorted_depths == unique_depths(u));
                x_sorted_depths_map=x_sorted_depths(idx_depths) ;
                y_sorted_depths_map=y_sorted_depths(idx_depths);
                z_sorted_depths_map=z_sorted_depths(idx_depths);
                non_adiabT_map{u}=depths_nonadiabatic_temperature_sorted(idx_depths);
                x_depths_downsampled{u} = x_sorted_depths_map(1:downsample_factor:end);
                y_depths_downsampled{u} = y_sorted_depths_map(1:downsample_factor:end);
                z_depths_downsampled{u} = z_sorted_depths_map(1:downsample_factor:end);
                if ~strcmp(additional_fields_to_load{1},'')
                    for f = 1:numel(additional_fields_to_load_sorted)
                        depths_additional_fields_to_load_map{u,f} = additional_fields_to_load_sorted{f}(idx_depths);
                    end
                end
            end
        end

        %Lets set griddata separately so we only call it only if needed
        if  any(strcmp(additional_postprocesses, 'subduction_and_plume_statistics')) || any(strcmp(output_additional_maps_figures, 'topography')) ...
                || strcmp(remove_subductions_to_oceanic_age, 'true')

            min_surface_topography = min(surface_topography_surface_downsampled);
            max_surface_topography = max(surface_topography_surface_downsampled);
            mean_surface_topography(count_ite) = mean(surface_topography_surface_downsampled);

            %         if init_step~=0
            %         %        Correction for mean topography subsiding
            %         if(count_ite<=Reference_mean_correction_steps)
            %             diff_mean_surface_topography(count_ite)=0;
            %             %             disp('test1')
            %         elseif (mean_surface_topography(count_ite)<mean_surface_topography(count_ite-1))
            %             diff_mean_surface_topography(count_ite)=mean_surface_topography(Reference_mean_correction_steps)-mean_surface_topography(count_ite);
            %             %             disp('test2')
            %         else
            %             diff_mean_surface_topography(count_ite)=0;
            %             %             disp('test3')
            %         end
            %         else
            % another way of correcting the topography in case we are not
            % starting from timestep 0
            if(count_ite_with_steps<=Reference_mean_correction_steps)
                diff_mean_surface_topography(count_ite)=0;
                %             disp('test1')
                %         elseif (mean_surface_topography(count_ite)<mean_surface_topography(count_ite-1))
            elseif (mean_surface_init_topography>mean_surface_topography(count_ite))
                diff_mean_surface_topography(count_ite)=mean_surface_init_topography-mean_surface_topography(count_ite);
                %             disp('test2')
            else
                diff_mean_surface_topography(count_ite)=0;
                %             disp('test3')
            end
            %         end

            % Reshape the surface topography data to match the grid size
            Topo = griddata(longitude_surface_downsampled, latitude_surface_downsampled, surface_topography_surface_downsampled, Xeq, Yeq);
            Topo = Topo(1:181,1:361)+diff_mean_surface_topography(count_ite);
        end


        % Check if 'Geofeatures' is called with 'subduction_and_plume_statistics'
        if strcmp(output_additional_maps_figures, 'geofeatures')
            if ~any(strcmp(additional_postprocesses, 'subduction_and_plume_statistics'))
                disp('Warning: "Geofeatures" mapping can only be call if "subduction_and_plume_statistics" is turned on');
            end
        end

        if  any(strcmp(additional_postprocesses, 'subduction_and_plume_statistics'))
            for uu = 1:numel(unique_depths)
                longitude_depths_downsampled = atan2(y_depths_downsampled{uu}, x_depths_downsampled{uu}) * 180 / pi;
                latitude_depths_downsampled = asin(z_depths_downsampled{uu} ./ sqrt(x_depths_downsampled{uu}.^2 + y_depths_downsampled{uu}.^2 + z_depths_downsampled{uu}.^2)) * 180 / pi;
                non_adiabT_depths{uu}=griddata(longitude_depths_downsampled, latitude_depths_downsampled, non_adiabT_map{uu}, Xeq, Yeq);
            end
        end

        if any(strcmp(additional_postprocesses, 'subduction_and_plume_statistics'))&& any(strcmp(output_additional_maps_figures, 'geofeatures'))||any(strcmp(output_additional_maps_figures, 'strain_rate'))
            strain_rate = griddata(longitude_surface_downsampled, latitude_surface_downsampled, surface_strain_rate_sorted, Xeq, Yeq);
        end

        if any(strcmp(files_requirement, 'lithosphere'))
            x_lithosphere = lithosphere.("Points:0");
            y_lithosphere = lithosphere.("Points:1");
            z_lithosphere = lithosphere.("Points:2");

            theta_lithosphere = atan2(y_lithosphere, x_lithosphere);
            [theta_sorted_lithosphere, theta_index_lithosphere] = sort(theta_lithosphere);
            x_sorted_lithosphere = x_lithosphere(theta_index_lithosphere);
            y_sorted_lithosphere = y_lithosphere(theta_index_lithosphere);
            z_sorted_lithosphere = z_lithosphere(theta_index_lithosphere);
            x_sorted_lithosphere = x_lithosphere(theta_index_lithosphere);
            y_sorted_lithosphere = y_lithosphere(theta_index_lithosphere);
            z_sorted_lithosphere = z_lithosphere(theta_index_lithosphere);
            oceanic_age_lithosphere_sorted = lithosphere.oceanic_age(theta_index_lithosphere);
            divergence_lithosphere_sorted = lithosphere.Divergence(theta_index_lithosphere);

            % load addtional variable to plot if required by the user
            if ~strcmp(additional_fields_to_load{1},'')
                additional_fields_to_load_map = cell(numel(additional_fields_to_load));
                for f = 1:numel(additional_fields_to_load)
                    if ismember(additional_fields_to_load{f}, lithosphere.Properties.VariableNames)
                        % Retrieve the data using a temporary variable
                        temp_data = lithosphere.(additional_fields_to_load{f})(theta_index_lithosphere);
                        % Assign to the cell array using curly braces
                        lithosphere_additional_fields_to_load_sorted{f} = temp_data;
                    end
                end
            end

            % Check if the field exists and assign the parameter accordingly
            if ismember(compositional_field_name_of_continents, lithosphere.Properties.VariableNames)
                lithosphere_continent_sorted = lithosphere.(compositional_field_name_of_continents)(theta_index_lithosphere);
            else
                % Handle the case where the field does not exist
                lithosphere_continent_sorted = 'NaN';  % or assign a default value
                % Display a message
                disp('Warning: There is no field assigned to continents. This is fine if your model has no continents, otherwise please assign a field so they can be tracked and projected.');
            end

            x_lithosphere_downsampled = x_sorted_lithosphere(1:downsample_factor:end);
            y_lithosphere_downsampled = y_sorted_lithosphere(1:downsample_factor:end);
            z_lithosphere_downsampled = z_sorted_lithosphere(1:downsample_factor:end);

            longitude_lithosphere_downsampled = atan2(y_lithosphere_downsampled, x_lithosphere_downsampled) * 180 / pi;
            latitude_lithosphere_downsampled = asin(z_lithosphere_downsampled ./ sqrt(x_lithosphere_downsampled.^2 + y_lithosphere_downsampled.^2 + z_lithosphere_downsampled.^2)) * 180 / pi;
            if any(strcmp(output_additional_maps_figures, 'divergence'))
                divergence = griddata(longitude_lithosphere_downsampled, latitude_lithosphere_downsampled, divergence_lithosphere_sorted, Xeq, Yeq);
            end
        end

        %Lets set griddata separately so we only call it only if needed
        if any(strcmp(additional_postprocesses, 'oceanic_age_statistics')) || any(strcmp(output_additional_maps_figures, 'oceanic_age'))
            oceanic_age = griddata(longitude_lithosphere_downsampled, latitude_lithosphere_downsampled, oceanic_age_lithosphere_sorted, Xeq, Yeq);
        end

        if any(strcmp(files_requirement, 'surface'))
            if ~strcmp(surface_continent_sorted,'NaN')
                continents = griddata(longitude_surface_downsampled, latitude_surface_downsampled, surface_continent_sorted, Xeq, Yeq);
            else
                %assign continent is null everywhere
                continents = zeros(size(Xeq));
            end
        elseif  ~any(strcmp(files_requirement, 'surface'))
            if any(strcmp(files_requirement, 'lithosphere'))
                disp('here2');
                if ~strcmp(surface_continent_sorted,'NaN')
                    continents = griddata(longitude_lithosphere_downsampled, latitude_lithosphere_downsampled, lithosphere_continent_sorted, Xeq, Yeq);
                else
                    %assign continent is null everywhere
                    continents = zeros(size(Xeq));
                end
            end
        end


        %if any additionnal field needs to be resampled
        if ~strcmp(additional_fields_to_load{1},'')
            disp('Resampling additional_fields_to_load')
            if any(strcmp(files_requirement, 'surface'))
                for f = 1:numel(additional_fields_to_load)
                    surface_additional_fields{f}= griddata(longitude_surface_downsampled, latitude_surface_downsampled, surface_additional_fields_to_load_sorted{f}, Xeq, Yeq);
                end
            end
            if  any(strcmp(files_requirement, 'lithosphere'))
                for f = 1:numel(additional_fields_to_load)
                    lithopshere_additional_fields = griddata(longitude_lithosphere_downsampled, latitude_lithosphere_downsampled, lithosphere_additional_fields_to_load_sorted{f}, Xeq, Yeq);
                end
            end
            if any(strcmp(files_requirement, 'depths'))
                for f = 1:numel(additional_fields_to_load)
                    for uu = 1:numel(unique_depths)
                        longitude_depths_downsampled = atan2(y_depths_downsampled{uu}, x_depths_downsampled{uu}) * 180 / pi;
                        latitude_depths_downsampled = asin(z_depths_downsampled{uu} ./ sqrt(x_depths_downsampled{uu}.^2 + y_depths_downsampled{uu}.^2 + z_depths_downsampled{uu}.^2)) * 180 / pi;
                        if ~isempty(additional_fields_depths_to_visualize_ite{f,uu})
                            depths_additional_fields{uu,f}=griddata(longitude_depths_downsampled, latitude_depths_downsampled, depths_additional_fields_to_load_map{uu,f}, Xeq, Yeq);
                        else
                            depths_additional_fields{uu,f}=[];
                        end
                    end
                end
            end
        end

        %% Statistics place
        if any(strcmp(additional_postprocesses, 'subduction_and_plume_statistics'))|| any(strcmp(output_additional_maps_figures, 'topography')) || any(strcmp(output_additional_maps_figures, 'geofeatures')) ...
                || any(strcmp(output_additional_maps_figures, 'strain_rate')) || any(strcmp(output_additional_maps_figures, 'oceanic_age')) || any(strcmp(additional_postprocesses, 'oceanic_age_statistics'))||~strcmp(additional_fields_to_load{1}, '')
            %% Plot topography corrected

            %continents positions
            continents_threshold = continents>0.3;
            idx_is_nan_to_zero=isnan(continents_threshold);
            continents_threshold(idx_is_nan_to_zero)=0;
            oceans = 1 - continents_threshold;
            %with resolution of Xeq and Yeq of 1deg from griddata
            idx_oceans=find(oceans == 1);
            idx_cont=find(oceans == 0);
            %         if i==1
            %             contour_continents_init=continents;
            %         end
            continents_reversed_position_for_plot = ind2rgb(oceans + 1, CMap_continents);
            continents_reversed_position_for_plot_dark = ind2rgb(oceans + 1, CMap_continents_dark);

        end

        if any(strcmp(additional_postprocesses, 'subduction_and_plume_statistics'))|| any(strcmp(output_additional_maps_figures, 'topography')) || any(strcmp(output_additional_maps_figures, 'geofeatures')) || strcmp(remove_subductions_to_oceanic_age, 'true')
            % Correct the topography for the water load in oceans
            smooth_Topo = imgaussfilt(Topo,0.3);
            Topo_corrected=smooth_Topo+topography_correction; %apply correction from bathymetry to sea levl
            Topo_corrected(idx_oceans) = Topo_corrected(idx_oceans) - (Topo_corrected(idx_oceans)<0).* (1000 / 3300).* abs(Topo_corrected(idx_oceans));
            %         abs(nanmean(Topo_corrected(idx_oceans)))+abs(nanmean(Topo_corrected(idx_cont)))

            %With 1900m readjustement due to wrong continent that needs some
            %adjustement
%             Topo_corrected(idx_cont) = Topo_corrected(idx_cont)+1900;
        end
        if any(strcmp(output_additional_maps_figures, 'topography'))
            colormap_name = 'Oleron_256';
            oleron_256 = load_modified_colors(colormap_name);
            num_lines  = size(oleron_256,1);

            Topo_color = interp1(oleron_256(:, 1),oleron_256(:, 2:end), linspace(min(oleron_256(:, 1)), max(oleron_256(:, 1)), num_lines));

            h=figure;
            % Set visibility based on Display_figures
            if strcmp(Display_figures, 'true')
                set(h, 'Visible', 'on');
            else
                set(h, 'Visible', 'off');
            end
            %         figure;
            set(gcf, 'Units', 'Normalized', 'OuterPosition', [0, 0.04, 1, 0.96],'color','w');clf;
            view(0, 90);
            axesm mercator;
            framem;
            gridm;
            ax0 = gca;
            setm(ax0, 'MapProjection', 'robinson');
            geoimg = geoshow(Yeq, Xeq,Topo_corrected,'DisplayType', 'texturemap'); % save object
            shadem(-14.5,[210 75])
            geoimg.AlphaDataMapping = 'none'; % interpet alpha values as transparency values
            geoimg.FaceAlpha = 'texturemap'; % Indicate that the transparency can be different each pixel
            alpha(geoimg,double(~isnan(Topo_corrected)))

            % Number of lines in the gradient
            colormap(ax0,Topo_color./num_lines); clim([-8500 4500]);
            c0=colorbar('Position',[0.1 0.18 0.02 0.66]);
            c0.Label.String= "Elevations [m]";
            c0.Label.FontSize = 14;
            hold on;
            ax2 = gca;
            setm(ax2, 'MapProjection', 'robinson');
            geoshow(Yeq, Xeq,continents_reversed_position_for_plot_dark,'facealpha', 0);
            contourm(Yeq,Xeq, continents, 0.3, 'k','LineWidth',2);
            contourm(Yeq,Xeq, contour_continents_init, 0.3, 'w','LineWidth',2);
            %         contourm(Yeq,Xeq, continents_init, 0.3, 'w','LineWidth',2);
            %         quivermc(Yeq_lr, Xeq_lr, V_theta_vec_surface, V_phi_vec_surface,'color','k','linewidth',1,'units','Velocity (m.yr^-^1)','reference',scaling_velocity_factor_surface);
            %             plume_str = ['Plumes = ' num2str(numPlumes)];
            %             subduction_str = ['Subductions = ' num2str(numSubductions)];
            %             title([{time_str} plume_str subduction_str]);
            title({time_str});
            if strcmp(plot_continents_border_from_reconstruction,'true')

                ax3 = gca;
                setm(ax3, 'MapProjection', 'robinson')
                %             geoshow(shape_continents, 'DisplayType', 'polygon', 'FaceColor', [0.61 0.38 0.32], 'EdgeColor', [0.41 0.38 0.32], 'LineWidth', 2, 'FaceAlpha', 0.1);
                geoshow(shape_continents, 'DisplayType', 'polygon', 'FaceColor', [0.61 0.38 0.32],'EdgeColor', 'none', 'LineWidth', 1, 'FaceAlpha', 0.15);

                ax3.SortMethod = 'childorder';
            end
            ax2.SortMethod='childorder';
            hold off;


            % Create the full path for the repository
            topography_output = fullfile(path_model_output, 'Maps_topography');

            % Check if the repository already exists
            if exist(topography_output, 'dir') ~= 7
                % Create the repository if it doesn't exist
                mkdir(topography_output);
                fprintf('Repository Maps_topography created successfully.');
            end

            % Save the figure in the repository
            fig_filename = fullfile(topography_output, sprintf('topography_%04d.png', i));
            saveas(gcf, fig_filename);
            fprintf('Figure saved: %s\n', fig_filename);

        end

        if any(strcmp(additional_postprocesses, 'subduction_and_plume_statistics')) && i<=min([num_files_surface, num_files_depths])
            %% Track plumes and subductions over depths
            for tt= 1:numel(plumes_depths_tracking)
                matching_depths_plumes_indices = find(ismember(unique_depths, plumes_depths_tracking(tt).*1e3));
                current_depth_plumes(tt)=unique_depths(matching_depths_plumes_indices);
                plumes_position = non_adiabT_depths{matching_depths_plumes_indices}>plumes_non_adiabatic_tracking_temperature(tt);
                plumes_positions2=1-plumes_position;
                %             boundaries_position3 = ind2rgb(Plate_boundaries2+1, CMap_boundaries);]
                cmatrix_plume = contourm(Yeq,Xeq, non_adiabT_depths{matching_depths_plumes_indices}, [plumes_non_adiabatic_tracking_temperature(tt) 99999], 'r','LineWidth',2);
                if (subduction_depths_tracking(tt)==0)
                    current_depth_subductions(tt) = 0;
                    subductions_position = Topo_corrected<=trenches_elevation_threshold;
                    trenches_position2=1-subductions_position;
                    cmatrix_subduction = contourm(Yeq,Xeq, Topo_corrected, [-90000 trenches_elevation_threshold], 'b','LineWidth',2);

                    %for plot
                    %                 subductions_position3 = ind2rgb(trenches_position2 + 1, CMap_subduction);
                elseif (subduction_depths_tracking(tt)~=0)
                    matching_depths_subductions_indices = find(ismember(unique_depths, subduction_depths_tracking(tt).*1e3));
                    current_depth_subductions(tt)=unique_depths(matching_depths_subductions_indices);
                    subductions_position = non_adiabT_depths{matching_depths_subductions_indices}<subduction_non_adiabatic_tracking_temperature(tt);
                    subductions_position2=1-subductions_position;
                    cmatrix_subduction = contourm(Yeq,Xeq, non_adiabT_depths{matching_depths_subductions_indices}, [-90000 subduction_non_adiabatic_tracking_temperature(tt)], 'b','LineWidth',2);

                    %For plot
                    %                 subductions_position3 = ind2rgb(subductions_position2 + 1, CMap_subduction);

                end
                %            cmatrix_boundaries = contourm(Yeq,Xeq, log10(strain_rate), [5e-15 1], 'm','LineWidth',2);
                %%
                if isempty(cmatrix_plume)
                    numPlumes(tt) = 0;
                else
                    [x_plume, y_plume] = C2xyz(cmatrix_plume);
                    polyin_plumes = polyshape(x_plume, y_plume);
                    numPlumes(tt) = polyin_plumes.NumRegions;
                end
                if isempty(cmatrix_subduction)
                    numSubductions(tt) = 0;
                    numSubduction_connected(tt) = 0;
                    Maxferet_mean(tt) = 0;
                    Maxferet_max(tt)= 0;
                    Maxferet_min(tt) = 0;
                    nferet(tt) = 0;
                else
                    [x_subduction, y_subduction] = C2xyz(cmatrix_subduction);
                    polyin_subductions = polyshape(x_subduction, y_subduction);
                    numSubductions(tt) = polyin_subductions.NumRegions;
                    It2 = bwmorph(subductions_position,'thin','inf');
                    %                 If some cleaning is needed
                    %                 It2 = bwmorph(It,'clean');
                    Subduction_connectivity= bwconncomp(It2);
                    numSubduction_connected(tt) = Subduction_connectivity.NumObjects;
                    %           Lets calculate the max "Feret" length of each subduction system even after connection
                    [out,LM] = bwferet(It2,'MaxFeretProperties');

                    % Use WGS84 ellipsoid model
                    wgs84 = wgs84Ellipsoid;
                    % Calculate individual subduction distances
                    maxLabel = max(LM(:));

                    % Plot the maximum feretdistance for each subduction zone
                    %             figure();
                    %             hf = imshow(LM,[]);
                    %             axis = hf.Parent;
                    %             for labelvalues = 1:maxLabel
                    %                 xmin = [out.MaxCoordinates{labelvalues}(1,1) out.MaxCoordinates{labelvalues}(2,1)];
                    %                 ymin = [out.MaxCoordinates{labelvalues}(1,2) out.MaxCoordinates{labelvalues}(2,2)];
                    %                 imdistline(axis,xmin,ymin);
                    %                 distance_feret_km(labelvalues) = distance(out.MaxCoordinates{labelvalues}(1,1)-180.5,out.MaxCoordinates{labelvalues}(1,2)-90.5,out.MaxCoordinates{labelvalues}(2,1)-180.5,out.MaxCoordinates{labelvalues}(2,2)-90.5,wgs84);
                    %                 distance_feret_km(labelvalues) = distance_feret_km(labelvalues)./1e3;
                    %             end
                    %             title(axis,'Maximum Feret Diameter of Objects');
                    %             colorbar('Ticks',1:maxLabel)

                    % We calculate the maximum distance from each endpoint of the subduction branches and convert to km
                    for labelvalues = 1:maxLabel
                        distance_feret_km(labelvalues) = distance(out.MaxCoordinates{labelvalues}(1,1)-180.5,out.MaxCoordinates{labelvalues}(1,2)-90.5,out.MaxCoordinates{labelvalues}(2,1)-180.5,out.MaxCoordinates{labelvalues}(2,2)-90.5,wgs84);
                        distance_feret_km(labelvalues) = distance_feret_km(labelvalues)./1e3;
                    end
                    Maxferet_mean(tt) = mean(distance_feret_km);
                    Maxferet_max(tt) = max(distance_feret_km);
                    Maxferet_min(tt) = min(distance_feret_km);
                    nferet(tt) = numel(distance_feret_km);

                    % %             We can plot the maximum ferret distribution for each subduction zone using an histogram
                    %             figure();
                    %             histogram(distance_feret_km,30);xlabel('Distance Feret (km)');
                    %             ylabel('Frequency');
                    %             title('Histogram of Distance Feret');

                    %             Not so usefull for now but this is a way to calculate the
                    %             circularity of the subductions, better will be to calculate
                    %             the radius of curvature of each subduction "trenches"
                    %             propriety_s = regionprops(It2,'Circularity');
                    %             Circularity_subduction = cat(1,propriety_s.Circularity)
                    %             [population2,gof] = fit(distance_feret_km',Circularity_subduction,'poly2')
                    %             figure();plot(population2,distance_feret_km',Circularity_subduction,'o')


                end

                %If the user want to plot an additional field with
                %threshold on top of the geofeature map then we run the
                %following. The addtional field  will be read from the
                %additional field function.
                
                CMap_additional_fields_for_geofeatures = [0.8 0 0.8 ;  1 1 1];

                if ~isempty(plot_an_additional_field_on_top_of_geofeatures)
                    if ismember(plot_an_additional_field_on_top_of_geofeatures, additional_fields_to_load)
                        additional_fields_threshold_ite_to_use_for_geofeatures = additional_fields_threshold_ite';
                        additional_fields_depths_to_visualize_ite_to_use_on_top = additional_fields_depths_to_visualize_ite';
                        [field_is_member, idx_member] = ismember(plot_an_additional_field_on_top_of_geofeatures, additional_fields_to_load);
                        if ~isempty(str2num(additional_fields_depths_to_visualize_ite_to_use_on_top{tt,idx_member}))
                            if ~isempty(additional_fields_threshold_ite_to_use_for_geofeatures{tt,idx_member})
                                if contains(additional_fields_threshold_ite_to_use_for_geofeatures{tt,idx_member}, '<')
                                    %                                             disp('test2')
                                    % Check if the string is not empty
                                    current_map_thlds = depths_additional_fields{tt,idx_member} < str2double(extractAfter(additional_fields_threshold_ite_to_use_for_geofeatures{tt,idx_member}, '<'));
                                    thld_fields = str2double(extractAfter(additional_fields_threshold_ite_to_use_for_geofeatures{tt,idx_member}, '<'));
                                else
                                    % Check if the string is not empty
                                    current_map_thlds = depths_additional_fields{tt,idx_member} > str2double(extractAfter(additional_fields_threshold_ite_to_use_for_geofeatures{tt,idx_member}, '>'));
                                    thld_fields = str2double(extractAfter(additional_fields_threshold_ite_to_use_for_geofeatures{tt,idx_member}, '>'));
                                end
                                if isnan(thld_fields)
                                    thld_fields = default_threshold_fields;
                                end
                                current_map_position2=1-current_map_thlds;
                                current_map_position3 = ind2rgb(current_map_position2 + 1, CMap_additional_fields_for_geofeatures);
                                current_map_used = current_map_position3;
                            end
                        end
                    end
                else
                    disp(['The field ', plot_an_additional_field_on_top_of_geofeatures, ' should be added in additional field preferably to plot for all depths to be plotted on top of geofeatures']);
                end


                if strcmp(output_figures_for_spherical_additional_postprocess,'true') && any(strcmp(output_additional_maps_figures, 'geofeatures'))
                    Plate_boundaries = strain_rate;
                    Plate_boundaries = Plate_boundaries>5e-15;
                    Plate_boundaries2 = 1-Plate_boundaries;

                    h=figure;
                    % Set visibility based on Display_figures
                    if strcmp(Display_figures, 'true')
                        set(h, 'Visible', 'on');
                    else
                        set(h, 'Visible', 'off');
                    end
                    set(gcf, 'Units', 'Normalized', 'OuterPosition', [0, 0.04, 1, 0.96],'color','w');clf;
                    set(gcf, 'color', 'w');
                    view(0, 90);
                    ax1 = gca;
                    axesm mercator;
                    framem;
                    gridm;
                    plumes_position3 = ind2rgb(plumes_positions2 + 1, CMap_plume);
                    geoshow(Yeq, Xeq, plumes_position3);
                    setm(ax1, 'MapProjection', 'robinson');
                    freezeColors(ax1);
                    hold on;
                    if (subduction_depths_tracking(tt)==0)
                        subductions_position3 = ind2rgb(trenches_position2 + 1, CMap_subduction);
                        cmatrix_subduction = contourm(Yeq,Xeq, Topo_corrected, [-90000 trenches_elevation_threshold], 'b','LineWidth',2);
                    elseif (subduction_depths_tracking(tt)~=0)
                        subductions_position3 = ind2rgb(subductions_position2 + 1, CMap_subduction);
                        cmatrix_subduction = contourm(Yeq,Xeq, non_adiabT_depths{matching_depths_subductions_indices}, [-90000 subduction_non_adiabatic_tracking_temperature(tt)], 'b','LineWidth',2);
                    end
                    geoshow(Yeq, Xeq,subductions_position3,'FaceAlpha',0.5);
                    boundaries_position3 = ind2rgb(Plate_boundaries2+1, CMap_boundaries);
                    %                     We can visualize strain rate only for the surface layers
                    %                     if(subduction_depths_tracking(tt)==0 || subduction_depths_tracking(tt)==440)
                    geoshow(Yeq, Xeq, boundaries_position3,'FaceAlpha',0.4);
                    %                     end
                    cmatrix_plume = contourm(Yeq,Xeq, non_adiabT_depths{matching_depths_plumes_indices}, [plumes_non_adiabatic_tracking_temperature(tt) plumes_non_adiabatic_tracking_temperature(tt)], 'r','LineWidth',2);
                    %                     plume_str = sprintf('Plumes: %s', strjoin(arrayfun(@(dtracking) sprintf('%d (%d km)', numPlumes(dtracking), plumes_depths_tracking(dtracking)), 1:numel(numPlumes), 'UniformOutput', false), ' '));
                    plume_str = sprintf('Plumes: %s', strjoin(arrayfun(@(dtracking) sprintf('%d (%d km, %dK)', numPlumes(dtracking), plumes_depths_tracking(dtracking), plumes_non_adiabatic_tracking_temperature(dtracking)), 1:numel(numPlumes), 'UniformOutput', false), ' '));
                    subduction_str = sprintf('Subductions: %s', strjoin(arrayfun(@(dtracking) sprintf('%d (%d km, %dK)', numSubductions(dtracking), subduction_depths_tracking(dtracking), subduction_non_adiabatic_tracking_temperature(dtracking)), 1:numel(numSubductions), 'UniformOutput', false), ' '));
                    %                     subduction_str = sprintf('Subductions: %s', strjoin(arrayfun(@(dtracking) sprintf('%d (%d km)', numSubductions(dtracking), subduction_depths_tracking(dtracking)), 1:numel(numSubductions), 'UniformOutput', false), ' '));
                    current_depth_str = sprintf('Current depth = %d km', current_depth_subductions(tt)./1e3);
                    title_str = [{time_str} current_depth_str plume_str subduction_str];
                    title(title_str);
                    contourm(Yeq,Xeq, contour_continents_init, 0.3, 'k','LineWidth',2);
                    %             contourm(Yeq,Xeq, continents_init, 0.3, 'k','LineWidth',2);
                    ax2 = gca;
                    setm(ax2, 'MapProjection', 'robinson');
                    geoshow(Yeq, Xeq,continents_reversed_position_for_plot,'facealpha', 0.3);
                    colormap(ax2);

                    if strcmp(plot_continents_border_from_reconstruction,'true')

                        ax3 = gca;
                        setm(ax3, 'MapProjection', 'robinson')
                        %             geoshow(shape_continents, 'DisplayType', 'polygon', 'FaceColor', [0.61 0.38 0.32], 'EdgeColor', [0.41 0.38 0.32], 'LineWidth', 2, 'FaceAlpha', 0.1);
                        geoshow(shape_continents, 'DisplayType', 'polygon', 'FaceColor', [0.61 0.38 0.32],'EdgeColor', 'none', 'LineWidth', 1, 'FaceAlpha', 0.15);

                        ax3.SortMethod = 'childorder';
                    end

                    if ~isempty(plot_an_additional_field_on_top_of_geofeatures)

                        ax4 = gca;
                        setm(ax4, 'MapProjection', 'robinson')
                        %             geoshow(shape_continents, 'DisplayType', 'polygon', 'FaceColor', [0.61 0.38 0.32], 'EdgeColor', [0.41 0.38 0.32], 'LineWidth', 2, 'FaceAlpha', 0.1);
%                         geoshow(shape_continents, 'DisplayType', 'polygon', 'FaceColor', [0.61 0.38 0.32],'EdgeColor', 'none', 'LineWidth', 1, 'FaceAlpha', 0.15);
                     geoshow(Yeq, Xeq, current_map_used,'FaceAlpha',0.4);
                     cmatrix_fields = contourm(Yeq, Xeq, depths_additional_fields{tt,idx_member}, [thld_fields thld_fields], 'LineColor', [0.8 0 0.8], 'LineWidth', 2);
                        
                     ax4.SortMethod = 'childorder';
                    end

                    ax2.SortMethod = 'childorder';
                    hold off;

                    % Create the full path for the repository
                    geofeatures_output = fullfile(path_model_output, 'Maps_geofeatures');

                    % Check if the repository already exists
                    if exist(geofeatures_output, 'dir') ~= 7
                        % Create the repository if it doesn't exist
                        mkdir(geofeatures_output);
                        fprintf('Repository Maps_geofeatures created successfully.');
                    end

                    % Save the figure in the repository
                    fig_filename = fullfile(geofeatures_output, sprintf('geofeatures_depth_%04d_%04d.png',current_depth_subductions(tt)./1e3, i));
                    saveas(gcf, fig_filename);
                    fprintf('Figure saved: %s\n', fig_filename);

                end
            end
        end

        %% Plot additional field selected by the user

        thld_fields = [];
        current_map = [];
        current_map_position3 = [];
        index_cmap = [];

        if ~strcmp(additional_fields_to_load{1}, '')
            additional_fields_threshold_ite_to_use = additional_fields_threshold_ite';
            additional_fields_depths_to_visualize_ite_to_use = additional_fields_depths_to_visualize_ite';
            non_empty_maps = {surface_additional_fields, lithosphere_additional_fields, depths_additional_fields};
            non_empty_indices = ~cellfun('isempty', non_empty_maps);

            if any(non_empty_indices)
                % Iterate over the non-empty maps and make plots
                for tt = 1:numel(non_empty_maps)
                                RMS_melt_depths=zeros(size(non_empty_maps{tt},1),1);
                                RMS_melt_surface=zeros(size(non_empty_maps{tt},1),1);
                                RMS_melt_sublithosphere=zeros(size(non_empty_maps{tt},1),1);
                    if non_empty_indices(tt)
                        for ii=1:numel(additional_fields_to_load)
                            for jj = 1:size(non_empty_maps{tt},1)
                                if ~isempty(str2num(additional_fields_depths_to_visualize_ite_to_use{jj,ii}))
                                    if ~isempty(additional_fields_threshold_ite_to_use{jj,ii})
                                        if contains(additional_fields_threshold_ite_to_use{jj,ii}, '<')
                                            %                                             disp('test2')
                                            % Check if the string is not empty
                                            current_map_thlds = non_empty_maps{tt}{jj,ii} < str2double(extractAfter(additional_fields_threshold_ite_to_use{jj,ii}, '<'));
                                            thld_fields = str2double(extractAfter(additional_fields_threshold_ite_to_use{jj,ii}, '<'));
                                        else
                                            % Check if the string is not empty
                                            current_map_thlds = non_empty_maps{tt}{jj,ii} > str2double(extractAfter(additional_fields_threshold_ite_to_use{jj,ii}, '>'));
                                            thld_fields = str2double(extractAfter(additional_fields_threshold_ite_to_use{jj,ii}, '>'));
                                        end
                                        if isnan(thld_fields)
                                            thld_fields = default_threshold_fields;
                                        end
                                        current_map_position2=1-current_map_thlds;
                                        current_map_position3 = ind2rgb(current_map_position2 + 1, CMap_additional_fields);
                                        index_cmap = 3;
                                        current_map_used = current_map_position3;

                                    else
                                        %                                         disp('test3')
                                        thld_fields = default_threshold_fields; %give a default value for contourm function
                                        current_map_field = non_empty_maps{tt}{jj,ii};
                                        index_cmap = 4;
                                        current_map_used = current_map_field;
                                    end

                                    h=figure;
                                    % Set visibility based on Display_figures
                                    if strcmp(Display_figures, 'true')
                                        set(h, 'Visible', 'on');
                                    else
                                        set(h, 'Visible', 'off');
                                    end
                                    set(gcf, 'Units', 'Normalized', 'OuterPosition', [0, 0.04, 1, 0.96]);clf;
                                    set(gcf, 'color', 'w');
                                    view(0, 90);
                                    axesm mercator;
                                    framem;
                                    gridm;

                                    ax1 = gca;
                                    setm(ax1, 'MapProjection', 'robinson');
                                    zlabel(additional_fields_to_load{ii});

                                    if index_cmap ==4
                                        if non_empty_indices(1)==tt %surface
                                            title_str = sprintf('%s (at surface)', additional_fields_to_load{ii});
                                        else
                                            title_str = sprintf('%s (%d km)', additional_fields_to_load{ii},unique_depths(1)./1e3);
                                        end
                                        title_str=strrep(title_str, '_', ' ');
                                        title({title_str, time_str});

                                        cmap1=crameri('lajolla');
                                        c1 = colorbar();%ax1,'Position',[0.91 0.18 0.02 0.66]
                                        try
                                            clim([min(min(current_map_used)) max(max(current_map_used))]);
                                        end
                                        %                                     clim([0 1]);
                                        set(ax1,'ColorScale');
                                        c1.Label.String = additional_fields_to_load{ii};
                                        c1.Label.FontSize = 14;
                                        colormap(ax1,cmap1);
                                        geoshow(Yeq, Xeq, current_map_used,'DisplayType', 'surface');
                                    elseif index_cmap ==3
                                        if non_empty_indices(1)==tt %surface
                                            title_str = sprintf('%s (%s at surface)', additional_fields_to_load{ii},num2str(thld_fields,'%.2f'));
                                        else
                                            title_str = sprintf('%s (%s at %d km)', additional_fields_to_load{ii},num2str(thld_fields,'%.2f'),unique_depths(jj)./1e3);
                                        end
                                        title_str=strrep(title_str, '_', ' ');
                                        title({title_str, time_str});

                                        geoshow(Yeq, Xeq, current_map_used,'FaceAlpha',0.4);
                                    end
                                    freezeColors(ax1);
                                    hold on;
                                    cmatrix_fields = contourm(Yeq, Xeq, non_empty_maps{tt}{jj,ii}, [thld_fields thld_fields], 'r', 'LineWidth', 2);
                                    contourm(Yeq,Xeq, contour_continents_init, 0.3, 'k','LineWidth',2);
                                    %             contourm(Yeq,Xeq, continents_init, 0.3, 'k','LineWidth',2);
                                    ax2 = gca;
                                    setm(ax2, 'MapProjection', 'robinson');
                                    geoshow(Yeq, Xeq,continents_reversed_position_for_plot,'facealpha', 0.3);
                                    colormap(ax2);

                                    if strcmp(plot_continents_border_from_reconstruction,'true')
                                        ax3 = gca;
                                        setm(ax3, 'MapProjection', 'robinson')
                                        %             geoshow(shape_continents, 'DisplayType', 'polygon', 'FaceColor', [0.61 0.38 0.32], 'EdgeColor', [0.41 0.38 0.32], 'LineWidth', 2, 'FaceAlpha', 0.1);
                                        geoshow(shape_continents, 'DisplayType', 'polygon', 'FaceColor', [0.61 0.38 0.32],'EdgeColor', 'none', 'LineWidth', 1, 'FaceAlpha', 0.15);

                                        ax3.SortMethod = 'childorder';
                                    end

                                    ax2.SortMethod = 'childorder';
                                    hold off;

                                    % Create the full path for the repository
                                    geofeatures_output = fullfile(path_model_output, 'additional_field_maps');

                                    % Check if the repository already exists
                                    if exist(geofeatures_output, 'dir') ~= 7
                                        % Create the repository if it doesn't exist
                                        mkdir(geofeatures_output);
                                        fprintf('Repository Maps_additional_fields created successfully.');
                                    end

                                    %                                 Save the figure in the repository
                                    if non_empty_indices(1)==tt %surface
                                        fig_filename = fullfile(geofeatures_output, sprintf('%s_depth_surface_%04d_threshold_%s_.png',additional_fields_to_load{ii}, i,num2str(thld_fields,'%.2f')));
                                    else
                                        fig_filename = fullfile(geofeatures_output, sprintf('%s_depth_%04d_%04d_threshold_%s_.png',additional_fields_to_load{ii},unique_depths(jj), i,num2str(thld_fields,'%.2f')));
                                    end
                                    saveas(gcf, fig_filename);
                                    fprintf('Figure saved: %s\n', fig_filename);

                                    if any(strcmp(additional_postprocesses, 'melt_statistics'))
                                        if strcmp(additional_fields_to_load{ii},'melt_fraction')
                                            if non_empty_indices(1)==tt %surface
                                               RMS_melt_surface(jj) = rms(non_empty_maps{tt}{jj,ii}(:),'omitnan');
                                            elseif non_empty_indices(2)==tt %surface
                                               RMS_melt_sublithosphere(jj) = rms(non_empty_maps{tt}{jj,ii}(:),'omitnan');
                                            elseif non_empty_indices(3)==tt %surface
                                               RMS_melt_depths(jj) = rms(non_empty_maps{tt}{jj,ii}(:),'omitnan');
                                            end
                                        end
                                    end
                                end
                            end
                        end
                    end
                end
            end
        end


        %%To remove once the code above works fine
        %         if ~strcmp(additional_fields_to_load{1}, '')
        %             non_empty_maps = {surface_additional_fields, lithosphere_additional_fields, depths_additional_fields};
        %             non_empty_indices = ~cellfun('isempty', non_empty_maps);
        %
        %             if any(non_empty_indices)
        %                 % Iterate over the non-empty maps and make plots
        %                 for tt = 1:numel(non_empty_maps)
        %                     if non_empty_indices(tt)
        %                         for ii=1:numel(additional_fields_to_load)
        %                             if (size(non_empty_maps{tt},1)==1)
        %                                 jj=1;
        %                                 if ~isempty(additional_fields_threshold{1})
        %                                     if contains(additional_fields_threshold{1}, '<')
        %                                         % Check if the string is not empty
        %                                         current_map_thlds = non_empty_maps{tt}{ii} < str2double(extractAfter(additional_fields_threshold{1}, '<'));
        %                                         thld_fields = str2double(extractAfter(additional_fields_threshold{1}, '<'));
        %                                     else
        %                                         % Check if the string is not empty
        %                                         current_map_thlds = non_empty_maps{tt}{ii} > str2double(extractAfter(additional_fields_threshold{1}, '>'));
        %                                         thld_fields = str2double(extractAfter(additional_fields_threshold{1}, '>'));
        %                                     end
        %                                     if isnan(thld_fields)
        %                                         thld_fields = default_threshold_fields;
        %                                     end
        %                                     current_map_position2=1-current_map_thlds;
        %                                     %                                     figure();pcolor(current_map_position3);shading flat;
        %
        %                                     %                                 1-current_map;
        %                                     current_map_position3 = ind2rgb(current_map_position2 + 1, CMap_additional_fields);
        %                                     index_cmap = 1;
        %                                     current_map_used = current_map_position3;
        %
        %                                 else
        %                                     thld_fields = default_threshold_fields; %give a default value for contourm function
        %                                     current_map = non_empty_maps{tt}{ii};
        %                                     index_cmap = 2;
        %                                     current_map_used = current_map;
        %                                 end
        %                                 %                                 current_map_used = current_map_position3;
        %                                 %                                 ~isempty(additional_fields_threshold{1}) * current_map_field + isempty(additional_fields_threshold{1}) * current_map_position3;
        %                                 h=figure;
%                                         % Set visibility based on Display_figures
%                                         if strcmp(Display_figures, 'true')
%                                             set(h, 'Visible', 'on');
%                                         else
%                                             set(h, 'Visible', 'off');
%                                         end
        %                                 set(gcf, 'Units', 'Normalized', 'OuterPosition', [0, 0.04, 1, 0.96]);clf;
        %                                 set(gcf, 'color', 'w');
        %                                 view(0, 90);
        %                                 axesm mercator
        %                                 framem;
        %                                 gridm;
        %
        %                                 ax1 = gca;
        %                                 setm(ax1, 'MapProjection', 'robinson')
        %                                 zlabel(additional_fields_to_load{ii});
        %                                 if index_cmap ==2
        %                                     title_str = sprintf('%s (%d km)', additional_fields_to_load{ii},unique_depths(1)./1e3);
        %                                     title({title_str, time_str});
        %                                     cmap1=crameri('lajolla');
        %                                     c1 = colorbar();%ax1,'Position',[0.91 0.18 0.02 0.66]
        %                                     clim([min(min(current_map_used)) max(max(current_map_used))]);
        %                                     %                                     clim([0 1]);
        %                                     set(ax1,'ColorScale');
        %                                     c1.Label.String = additional_fields_to_load{ii};
        %                                     c1.Label.FontSize = 14;
        %                                     colormap(ax1,cmap1);
        %                                     geoshow(Yeq, Xeq, current_map_used,'DisplayType', 'surface');
        %
        %                                 else
        %                                     title_str = sprintf('%s (%s at %d km)', additional_fields_to_load{ii},num2str(thld_fields,'%.2f'),unique_depths(1)./1e3);
        %                                     title({title_str, time_str});
        %                                     geoshow(Yeq, Xeq, current_map_used,'FaceAlpha',0.4);
        %                                 end
        %
        %                                 freezeColors(ax1);
        %                                 hold on;
        % %                                 cmatrix_fields = contourm(Yeq, Xeq, non_empty_maps{tt}{ii}, [thld_fields thld_fields], 'r', 'LineWidth', 2);
        %                                 cmatrix_fields = contourm(Yeq, Xeq, non_empty_maps{tt}{ii},thld_fields, 'r', 'LineWidth', 2);
        %
        %                                 contourm(Yeq,Xeq, contour_continents_init, 0.3, 'k','LineWidth',2);
        %                                 %             contourm(Yeq,Xeq, continents_init, 0.3, 'k','LineWidth',2);
        %                                 ax2 = gca;
        %                                 setm(ax2, 'MapProjection', 'robinson')
        %                                 geoshow(Yeq, Xeq,continents_reversed_position_for_plot,'facealpha', 0.3);
        %                                 colormap(ax2);
        %
        %                                 if strcmp(plot_continents_border_from_reconstruction,'true')
        %
        %                                     ax3 = gca;
        %                                     setm(ax3, 'MapProjection', 'robinson')
        %                                     %             geoshow(shape_continents, 'DisplayType', 'polygon', 'FaceColor', [0.61 0.38 0.32], 'EdgeColor', [0.41 0.38 0.32], 'LineWidth', 2, 'FaceAlpha', 0.1);
        %                                     geoshow(shape_continents, 'DisplayType', 'polygon', 'FaceColor', [0.61 0.38 0.32],'EdgeColor', 'none', 'LineWidth', 1, 'FaceAlpha', 0.15);
        %
        %                                     ax3.SortMethod = 'childorder';
        %                                 end
        %
        %                                 ax2.SortMethod = 'childorder';
        %
        %                                 % Create the full path for the repository
        %                                 geofeatures_output = fullfile(path_model_output, 'additional_field_maps');
        %
        %                                 % Check if the repository already exists
        %                                 if exist(geofeatures_output, 'dir') ~= 7
        %                                     % Create the repository if it doesn't exist
        %                                     mkdir(geofeatures_output);
        %                                     fprintf('Repository Maps_additional_fields created successfully.');
        %                                 end
        %
        %                                 %                                 Save the figure in the repository
        %                                 fig_filename = fullfile(geofeatures_output, sprintf('%s_depth_%04d_%04d_threshold_%s_.png',additional_fields_to_load{ii},unique_depths(jj), i,num2str(thld_fields,'%.2f')));
        %                                 saveas(gcf, fig_filename);
        %                                 fprintf('Figure saved: %s\n', fig_filename);
        %                             else
        %                                 for jj = 1:size(non_empty_maps{tt},1)
        %                                     if ~isempty(additional_fields_threshold{tt})
        %                                         if contains(additional_fields_threshold{tt}, '<')
        %                                             disp('test22')
        %                                             % Check if the string is not empty
        %                                             current_map_thlds = non_empty_maps{tt}{jj,ii} < str2double(extractAfter(additional_fields_threshold{jj}, '<'));
        %                                             thld_fields = str2double(extractAfter(additional_fields_threshold{jj}, '<'));
        %                                         else
        %                                             % Check if the string is not empty
        %                                             current_map_thlds = non_empty_maps{tt}{jj,ii} > str2double(extractAfter(additional_fields_threshold{jj}, '>'));
        %                                             thld_fields = str2double(extractAfter(additional_fields_threshold{jj}, '>'));
        %                                         end
        %                                         if isnan(thld_fields)
        %                                             thld_fields = default_threshold_fields;
        %                                         end
        %                                         current_map_position2=1-current_map_thlds;
        %                                         current_map_position3 = ind2rgb(current_map_position2 + 1, CMap_additional_fields);
        %                                         index_cmap = 3;
        %                                         current_map_used = current_map_position3;
        %
        %                                     else
        %                                         disp('test22')
        %                                         thld_fields = default_threshold_fields; %give a default value for contourm function
        %                                         current_map_field = non_empty_maps{tt}{jj,ii};
        %                                         index_cmap = 4;
        %                                         current_map_used = current_map_field;
        %                                     end
        %
        %                                     h=figure;
%                                                 % Set visibility based on Display_figures
%                                                 if strcmp(Display_figures, 'true')
%                                                     set(h, 'Visible', 'on');
%                                                 else
%                                                     set(h, 'Visible', 'off');
%                                                 end
        %                                     set(gcf, 'Units', 'Normalized', 'OuterPosition', [0, 0.04, 1, 0.96]);clf;
        %                                     set(gcf, 'color', 'w');
        %                                     view(0, 90);
        %                                     axesm mercator
        %                                     framem;
        %                                     gridm;
        %
        %                                     ax1 = gca;
        %                                     setm(ax1, 'MapProjection', 'robinson')
        %                                     zlabel(additional_fields_to_load{ii});
        %
        %                                     if index_cmap ==4
        %                                         title_str = sprintf('%s (%d km)', additional_fields_to_load{ii},unique_depths(1)./1e3);
        %                                         title({title_str, time_str});
        %                                         cmap1=crameri('lajolla');
        %                                         c1 = colorbar();%ax1,'Position',[0.91 0.18 0.02 0.66]
        %                                         clim([min(min(current_map_used)) max(max(current_map_used))]);
        %                                         %                                     clim([0 1]);
        %                                         set(ax1,'ColorScale');
        %                                         c1.Label.String = additional_fields_to_load{ii};
        %                                         c1.Label.FontSize = 14;
        %                                         colormap(ax1,cmap1);
        %                                         geoshow(Yeq, Xeq, current_map_used,'DisplayType', 'surface');
        %                                     elseif index_cmap ==3
        %                                         title_str = sprintf('%s (%s at %d km)', additional_fields_to_load{ii},num2str(thld_fields,'%.2f'),unique_depths(jj)./1e3);
        %                                         title({title_str, time_str});
        %                                         geoshow(Yeq, Xeq, current_map_used,'FaceAlpha',0.4);
        %                                     end
        %                                     freezeColors(ax1);
        %                                     hold on;
        %                                     cmatrix_fields = contourm(Yeq, Xeq, non_empty_maps{tt}{jj,ii}, [thld_fields thld_fields], 'r', 'LineWidth', 2);
        %                                     contourm(Yeq,Xeq, contour_continents_init, 0.3, 'k','LineWidth',2);
        %                                     %             contourm(Yeq,Xeq, continents_init, 0.3, 'k','LineWidth',2);
        %                                     ax2 = gca;
        %                                     setm(ax2, 'MapProjection', 'robinson')
        %                                     geoshow(Yeq, Xeq,continents_reversed_position_for_plot,'facealpha', 0.3);
        %                                     colormap(ax2);
        %
        %                                     if strcmp(plot_continents_border_from_reconstruction,'true')
        %                                         ax3 = gca;
        %                                         setm(ax3, 'MapProjection', 'robinson')
        %                                         %             geoshow(shape_continents, 'DisplayType', 'polygon', 'FaceColor', [0.61 0.38 0.32], 'EdgeColor', [0.41 0.38 0.32], 'LineWidth', 2, 'FaceAlpha', 0.1);
        %                                         geoshow(shape_continents, 'DisplayType', 'polygon', 'FaceColor', [0.61 0.38 0.32],'EdgeColor', 'none', 'LineWidth', 1, 'FaceAlpha', 0.15);
        %
        %                                         ax3.SortMethod = 'childorder';
        %                                     end
        %
        %                                     ax2.SortMethod = 'childorder';
        %
        %                                     % Create the full path for the repository
        %                                     geofeatures_output = fullfile(path_model_output, 'additional_field_maps');
        %
        %                                     % Check if the repository already exists
        %                                     if exist(geofeatures_output, 'dir') ~= 7
        %                                         % Create the repository if it doesn't exist
        %                                         mkdir(geofeatures_output);
        %                                         fprintf('Repository Maps_additional_fields created successfully.');
        %                                     end
        %
        %                                     %                                 Save the figure in the repository
        %                                     fig_filename = fullfile(geofeatures_output, sprintf('%s_depth_%04d_%04d_threshold_%s_.png',additional_fields_to_load{ii},unique_depths(jj), i,num2str(thld_fields,'%.2f')));
        %                                     saveas(gcf, fig_filename);
        %                                     fprintf('Figure saved: %s\n', fig_filename);
        %
        %                                 end
        %                             end
        %                         end
        %                     end
        %                 end
        %             end
        %         end



        if strcmp(output_figures_for_spherical_additional_postprocess,'true') && any(strcmp(output_additional_maps_figures, 'strain_rate'))
            %% Plot Strain rate
            h=figure;
            % Set visibility based on Display_figures
            if strcmp(Display_figures, 'true')
                set(h, 'Visible', 'on');
            else
                set(h, 'Visible', 'off');
            end
            set(gcf, 'Units', 'Normalized', 'OuterPosition', [0, 0.04, 1, 0.96]);clf;
            set(gcf, 'color', 'w');
            view(0, 90);
            axesm mercator;

            geoshow(Yeq, Xeq, log10(strain_rate), 'DisplayType', 'surface');
            ax1 = gca;
            setm(ax1, 'MapProjection', 'robinson');
            zlabel('Surface');
            title({'Surface', time_str});
            c = colorbar;
            cmap = crameri('-roma');
            colormap(ax1,cmap);clim(ax1,[-16 -14]);
            set(ax1,'ColorScale','log');
            ticks = linspace(log10(1e-16), log10(1e-14), 6);
            tickLabels = arrayfun(@(x) sprintf('{%0.1f}', x), ticks, 'UniformOutput', false);
            c.Ticks = ticks;
            c.TickLabels = tickLabels;
            c.Label.String = 'Strain rate [1e^ /s]';
            freezeColors;
            framem;
            gridm;
            hold on;
            % Plot continents with transparency
            ax2 = gca;
            setm(ax2, 'MapProjection', 'robinson');
            geoshow(Yeq, Xeq,continents_reversed_position_for_plot,'facealpha', 0.3);
            colormap(ax2);freezeColors;
            contourm(Yeq,Xeq, continents, 0.3, 'k','LineWidth',2);
            contourm(Yeq,Xeq, contour_continents_init, 0.3, 'w','LineWidth',2);

            %         quivermc(Yeq_lr, Xeq_lr, V_theta_vec_surface, V_phi_vec_surface,'color','k','linewidth',1,'units','Velocity (m.yr^-^1)','reference',scaling_velocity_factor_surface);
            if strcmp(plot_continents_border_from_reconstruction,'true')

                ax3 = gca;
                setm(ax3, 'MapProjection', 'robinson')
                %             geoshow(shape_continents, 'DisplayType', 'polygon', 'FaceColor', [0.61 0.38 0.32], 'EdgeColor', [0.41 0.38 0.32], 'LineWidth', 2, 'FaceAlpha', 0.1);
                geoshow(shape_continents, 'DisplayType', 'polygon', 'FaceColor', [1 1 1],'EdgeColor', 'none', 'LineWidth', 1, 'FaceAlpha', 0.5);

                ax3.SortMethod = 'childorder';
            end

            ax2.SortMethod='childorder';
            hold off;

            % Create the full path for the repository
            strain_rate_output = fullfile(path_model_output, 'Maps_strain_rate');

            % Check if the repository already exists
            if exist(strain_rate_output, 'dir') ~= 7
                % Create the repository if it doesn't exist
                mkdir(strain_rate_output);
                fprintf('Repository Maps_strain_rate created successfully.');
            end

            % Save the figure in the repository
            fig_filename = fullfile(strain_rate_output, sprintf('strain_rate_%04d.png', i));
            saveas(gcf, fig_filename);
            fprintf('Figure saved: %s\n', fig_filename);


        end

        if any(strcmp(additional_postprocesses, 'oceanic_age_statistics')) || any(strcmp(output_additional_maps_figures, 'oceanic_age'))
            %% Plot oceanic age
            %                 idx_oceanic_domain = find(oceanic == 0);
            oceanic_age_domain_only=oceanic_age;
            oceanic_age_domain_only(continents_threshold) = NaN;
            oceanic_age_domain =oceanic_age_domain_only;

            if strcmp(remove_subductions_to_oceanic_age, 'true')

                subduction_zones = find(Topo_corrected<trenches_elevation_threshold);
                oceanic_age_domain_only(subduction_zones) = NaN;
            end

        end
        %         Get the active plate boundaries using a threshold on strain rate

        %         Plot similar to Coltice et al 2012, distribution of oceanic age

        if any(strcmp(additional_postprocesses, 'oceanic_age_statistics'))
            bin_edges_oce_age = 0:1:400;

            % Histogram
            h=figure;
            % Set visibility based on Display_figures
            if strcmp(Display_figures, 'true')
                set(h, 'Visible', 'on');
            else
                set(h, 'Visible', 'off');
            end
            set(gcf, 'color', 'w');
            yyaxis left;
            bin_counts = histcounts(oceanic_age_domain_only, bin_edges_oce_age);
            histogram(oceanic_age_domain_only, bin_edges_oce_age, 'FaceColor', 'blue', 'EdgeColor', 'black');hold on;
            ylim([0, 1200]);

            % Area per unit of age
            yyaxis right;
            bin_centers = (bin_edges_oce_age(1:end-1) + bin_edges_oce_age(2:end)) / 2;
            area_per_unit_of_age = bin_counts .* 111 * 111;
            plot(bin_centers, area_per_unit_of_age./1e6, 'LineWidth', 2, 'Color', 'red');

            % Kernel smoothing fit
            kernel_density = fitdist(oceanic_age_domain_only(:), 'Kernel', 'Kernel', 'epanechnikov', 'Bandwidth', 0.3);

            x_values = linspace(0, 250, 250);
            %         m = mean(kernel_density)
            %         med = median(kernel_density)
            %         s = std(kernel_density)
            smoothed_oceanic_age_density = pdf(kernel_density, x_values);
            %         I gave the  smoothed_oceanic_age_densit *1e3 just for the plot, the original value will be sent to the stats
            plot(x_values, smoothed_oceanic_age_density.*1e3, 'LineWidth', 2, 'Color', 'green');

            xlabel('Oceanic Age (My)');
            yyaxis left;ylabel('Frequency');
            yyaxis right;ylabel('Area per unit of age [.1e6 km^2 yr^{-1}]');
            title(['Oceanic age distribution at ', {time_str}]);
            legend('Age frequency', 'Area per unit of age','Smoothed Probability density function');

            % Create the full path for the repository
            oceanic_age_output = fullfile(path_model_output, 'Maps_oceanic_age');

            % Check if the repository already exists
            if exist(oceanic_age_output, 'dir') ~= 7
                % Create the repository if it doesn't exist
                mkdir(oceanic_age_output);
                fprintf('Repository Maps_oceanic_age created successfully.');
            end

            % Save the figure in the repository
            fig_filename = fullfile(oceanic_age_output, sprintf('oceanic_age_distribution_%04d.png', i));
            saveas(gcf, fig_filename);
            fprintf('Figure saved: %s\n', fig_filename);



        end

        if strcmp(output_figures_for_spherical_additional_postprocess,'true') && any(strcmp(output_additional_maps_figures, 'oceanic_age'))
            h=figure;
            % Set visibility based on Display_figures
            if strcmp(Display_figures, 'true')
                set(h, 'Visible', 'on');
            else
                set(h, 'Visible', 'off');
            end
            set(gcf, 'Units', 'Normalized', 'OuterPosition', [0, 0.04, 1, 0.96]);clf;
            set(gcf, 'color', 'w');
            view(0, 90);
            axesm mercator;
            framem;
            gridm;

            geoshow(Yeq, Xeq, oceanic_age_domain, 'DisplayType', 'surface');
            ax1 = gca;
            setm(ax1, 'MapProjection', 'robinson');
            zlabel('Surface');
            title({'Surface', time_str});
            cmap1=crameri('-roma');
            colormap(ax1,cmap1);
            c1 = colorbar(ax1,'Position',[0.91 0.18 0.02 0.66]);
            clim([0 200]);
            set(ax1,'ColorScale');
            ticks = linspace(0, 400, 6);
            tickLabels = arrayfun(@(x) sprintf('{%0.1f}', x), ticks, 'UniformOutput', false);
            c1.Ticks = ticks;
            c1.TickLabels = tickLabels;
            c1.Label.String = 'Oceanice age [My]';
            c1.Label.FontSize = 14;
            freezeColors;


            % freezeColors(c1)
            % Plot continents with transparency
            ax2 = gca;
            setm(ax2, 'MapProjection', 'robinson');
            geoshow(Yeq, Xeq, continents_reversed_position_for_plot,'facealpha', 0.3);
            colormap(ax2);
            contourm(Yeq,Xeq, continents, 0.3, 'k','LineWidth',2);
            contourm(Yeq,Xeq, contour_continents_init, 0.3, 'w','LineWidth',2);

            if strcmp(remove_subductions_to_oceanic_age, 'true')
                subduction_zones = contourm(Yeq,Xeq, Topo_corrected, [-90000 trenches_elevation_threshold], 'b','LineWidth',2);
            end
            if strcmp(plot_continents_border_from_reconstruction,'true')

                ax3 = gca;
                setm(ax3, 'MapProjection', 'robinson')
                %             geoshow(shape_continents, 'DisplayType', 'polygon', 'FaceColor', [0.61 0.38 0.32], 'EdgeColor', [0.41 0.38 0.32], 'LineWidth', 2, 'FaceAlpha', 0.1);
                geoshow(shape_continents, 'DisplayType', 'polygon', 'FaceColor', [0.61 0.38 0.32],'EdgeColor', 'none', 'LineWidth', 1, 'FaceAlpha', 0.15);

                ax3.SortMethod = 'childorder';
            end

            %         quivermc(Yeq_lr, Xeq_lr, V_theta_vec_surface, V_phi_vec_surface,'color','k','linewidth',1,'units','Velocity (m.yr^-^1)','reference',scaling_velocity_factor_surface);
            ax2.SortMethod='childorder';
            hold off;


            % Create the full path for the repository
            oceanic_age_output = fullfile(path_model_output, 'Maps_oceanic_age');

            % Check if the repository already exists
            if exist(oceanic_age_output, 'dir') ~= 7
                % Create the repository if it doesn't exist
                mkdir(oceanic_age_output);
                fprintf('Repository Maps_oceanic_age created successfully.');
            end

            % Save the figure in the repository
            fig_filename = fullfile(oceanic_age_output, sprintf('oceanic_age_%04d.png', i));
            saveas(gcf, fig_filename);
            fprintf('Figure saved: %s\n', fig_filename);


        end

        %% Plot divergence
        if strcmp(output_figures_for_spherical_additional_postprocess,'true') && any(strcmp(output_additional_maps_figures, 'divergence'))

            %Smooth the noisy divergence without loosing too much of the signal
            divergence2 = divergence;
            smooth_divergence = imgaussfilt(divergence2, 1);



            %Potentially we can remove low values of divergence
            % Find the maximum and minimum absolute values
            %         max_abs_value = max(smooth_divergence(:));
            %         min_abs_value = min(smooth_divergence(:));
            %         mean_abs_value = nanmean(smooth_divergence(:));
            % Set different thresholds for positive and negative values
            % threshold_pos = max_abs_value * 0.10;
            % threshold_neg = min_abs_value * 0.02;
            %

            % Different scheme can be use to choose the intervale of values to remove
            %         threshold_pos = (max_abs_value-mean_abs_value/2)*0.25;
            %         % * 0.10;
            %         threshold_neg = (min_abs_value-mean_abs_value)/2*0.25 ;

            %         threshold_pos = mode_abs_value-mode_abs_value*0.25;
            %         % * 0.10;
            %         threshold_neg = mode_abs_value+mode_abs_value*0.25 ;
            %         % * 0.1;
            %
            %         % Apply threshold conditions and set values to NaN
            %         smooth_divergence(smooth_divergence < threshold_pos & smooth_divergence> threshold_neg) = NaN;
            %

            %This plot check over a section how good is the smoothing compared
            %to the real divergence that is too noisy
            %         figure;
            %
            %         % Plot original data
            %         plot(divergence(:, 180), 'LineWidth', 2, 'DisplayName', 'Original Divergence');
            %         hold on;
            %
            %         % Plot smoothed data
            %         plot(smooth_divergence(:, 180), 'LineWidth', 2, 'DisplayName', 'Smoothed Divergence');
            %
            %         % Customize plot
            %         xlabel('Index');
            %         ylabel('Divergence');
            %         title('Comparison of Original and Smoothed Divergence');
            %         legend('show');set(gcf, 'color', 'w');
            %         grid on;


            h=figure;
            % Set visibility based on Display_figures
            if strcmp(Display_figures, 'true')
                set(h, 'Visible', 'on');
            else
                set(h, 'Visible', 'off');
            end
            set(gcf, 'Units', 'Normalized', 'OuterPosition', [0, 0.04, 1, 0.96]);clf;
            set(gcf, 'color', 'w');
            view(0, 90);
            axesm mercator;
            framem;
            gridm;

            geoshow(Yeq, Xeq,smooth_divergence, 'DisplayType', 'surface');
            ax1 = gca;
            setm(ax1, 'MapProjection', 'robinson');
            zlabel('Surface');
            title({'Surface', time_str});
            cmap1=crameri('vik');
            colormap(ax1,cmap1)
            c1 = colorbar();%ax1,'Position',[0.91 0.18 0.02 0.66]
            clim([-1e-8 1e-8]);
            set(ax1,'ColorScale');
            % ticks = linspace((-5e-7), (5e-7), 6);
            % tickLabels = arrayfun(@(x) sprintf('{%0.1f}', x), ticks, 'UniformOutput', false);
            % c1.Ticks = ticks;
            % c1.TickLabels = tickLabels;
            c1.Label.String = 'Divergence [1e^ /s]';
            c1.Label.FontSize = 14;
            freezeColors;



            ax2 = gca;
            setm(ax2, 'MapProjection', 'robinson');
            geoshow(Yeq, Xeq, continents_reversed_position_for_plot,'facealpha', 0.3);
            colormap(ax2);
            contourm(Yeq,Xeq, continents, 0.3, 'k','LineWidth',2);
            contourm(Yeq,Xeq, contour_continents_init, 0.3, 'w','LineWidth',2)
            %         quivermc(Yeq_lr, Xeq_lr, V_theta_vec_surface, V_phi_vec_surface,'color','k','linewidth',1,'units','Velocity (m.yr^-^1)','reference',scaling_velocity_factor_surface);
            if strcmp(plot_continents_border_from_reconstruction,'true')

                ax3 = gca;
                setm(ax3, 'MapProjection', 'robinson')
                %             geoshow(shape_continents, 'DisplayType', 'polygon', 'FaceColor', [0.61 0.38 0.32], 'EdgeColor', [0.41 0.38 0.32], 'LineWidth', 2, 'FaceAlpha', 0.1);
                geoshow(shape_continents, 'DisplayType', 'polygon', 'FaceColor', [0.61 0.38 0.32],'EdgeColor', 'none', 'LineWidth', 1, 'FaceAlpha', 0.15);

                ax3.SortMethod = 'childorder';
            end
            ax2.SortMethod='childorder';
            hold off;
            %%

            % Create the full path for the repository
            divergence_output = fullfile(path_model_output, 'Maps_divergence');

            % Check if the repository already exists
            if exist(divergence_output, 'dir') ~= 7
                % Create the repository if it doesn't exist
                mkdir(divergence_output);
                fprintf('Repository Maps_divergence created successfully.');
            end

            % Save the figure in the repository
            fig_filename = fullfile(divergence_output, sprintf('divergence_%04d.png', i));
            saveas(gcf, fig_filename);
            fprintf('Figure saved: %s\n', fig_filename);
        end



        %         if strcmp(output_figures_for_spherical_additional_postprocess,'true') && any(strcmp(output_additional_maps_figures, 'tectonics_regime'))
        %
        %
        %             % Assuming the relevant columns in your table are named principal_stress_direction_10, principal_stress_direction_11, and principal_stress_direction_12
        %
        %             % Extract the principal stress direction vectors
        %             S1 = [surface.principal_stress_direction_10, surface.principal_stress_direction_11, surface.principal_stress_direction_12];
        %
        %             % Compute the plunge angles for each vector
        %             plunge_angles = atan2d(sqrt(S1(:,1).^2 + S1(:,2).^2), abs(S1(:,3)));
        %
        %             % Define thresholds for classification
        %             normal_faulting_threshold = 30; % You may need to adjust this threshold based on your data
        %             strike_slip_threshold = 60; % You may need to adjust this threshold based on your data
        %
        %             % Initialize tectonic regime column
        %             surface.TectonicRegime = cell(size(surface, 1), 1);
        %
        %             % Classify based on plunge angles
        %             surface.TectonicRegime(plunge_angles <= normal_faulting_threshold) = {'NF'};
        %             surface.TectonicRegime(plunge_angles > normal_faulting_threshold & plunge_angles <= strike_slip_threshold) = {'NS'};
        %             surface.TectonicRegime(plunge_angles > strike_slip_threshold) = {'SS'};
        %
        %             % Display the result
        %             disp(surface(:, {'TectonicRegime'}));
        %         end


        %% Statistic for Vrms of continent
        if any(strcmp(additional_postprocesses, 'continents_VRMS')) && i<=num_files_surface
            %Lets avoid to use griddata, lets calculate the index
            idx_continents_sorted = find(surface_continent_sorted>0.3);
            norm_square_vel_surface = sqrt(surface_velocity0_sorted.^2+surface_velocity1_sorted.^2+surface_velocity2_sorted.^2);
            norm_square_vel_continents = sqrt(surface_velocity0_sorted(idx_continents_sorted).^2+surface_velocity1_sorted(idx_continents_sorted).^2+surface_velocity2_sorted(idx_continents_sorted).^2);
            V_rms_surface = rms(norm_square_vel_surface);
            V_rms_continents = rms(norm_square_vel_continents);
        end


        if strcmp(write_geofeatures_statistics, 'true')

            if any(strcmp(additional_postprocesses, 'subduction_and_plume_statistics'))
                ite_depth = numel(plumes_depths_tracking);
            else
                ite_depth = 1; 
            end
                if any(strcmp(additional_postprocesses, 'continents_VRMS')) || any(strcmp(additional_postprocesses, 'subduction_and_plume_statistics'))|| any(strcmp(additional_postprocesses, 'melt_statistics'))
                %% Write everything down
                % We already check in postprocess.m whether a file already exists, so that
                % we plot the data directly from the file and do not rewrite the file.
                % However, if someone wants to add data to the same file again
                % and again plotting can be switched off with postprocess.m
                % and the file could be opened and written again, thanks to the first
                % condition ~exist('geofeatures', 'var')

                filename_geofeatures = ['Geofeatures_' model_title '.csv'];
                fullFilePath_geofeatures = fullfile(path_model_output, filename_geofeatures);

                %     % Check if the file already exists
                %     if exist(fullFilePath_geofeatures, 'file')
                %         % File exists, load the existing table
                %         geofeatures = readtable(fullFilePath_geofeatures);
                %
                %         % Add new values to the existing table
                %         newRow = table(time_for_output_structures, numPlumes, numSubductions, numSubduction_connected, current_depth./1e3, Maxferet_mean, Maxferet_max, Maxferet_min, nferet, V_rms_surface, V_rms_continents, ...
                %             'VariableNames', {'Time', 'NumPlumes', 'NumSubductions', 'NumSubduction_connected', 'Depths', 'Maxferet_mean', 'Maxferet_max', 'Maxferet_in', 'nferet', 'V_rms_surface', 'V_rms_continents'});
                %
                %         geofeatures = vertcat(geofeatures, newRow);
                %     else
                %         % File doesn't exist, create a new table
                %         geofeatures = table(time_for_output_structures, numPlumes, numSubductions, numSubduction_connected, current_depth./1e3, Maxferet_mean, Maxferet_max, Maxferet_min, nferet, V_rms_surface, V_rms_continents,...
                %             'VariableNames', {'Time', 'NumPlumes', 'NumSubductions', 'NumSubduction_connected', 'Depths', 'Maxferet_mean', 'Maxferet_max', 'Maxferet_in', 'nferet', 'V_rms_surface', 'V_rms_continents'});
                %     end
                %
                %         %         % Display the updated table
                %     disp(geofeatures);
                %
                %     writetable(geofeatures, fullFilePath_geofeatures);
                %     fprintf('Table saved to %s successfully.\n', filename_geofeatures);
                %
                %We just want a depth
                current_depth = current_depth_subductions;


                if any(isnan(RMS_melt_surface)) || any(isnan(RMS_melt_sublithosphere)) || any(isnan(RMS_melt_depths))
                    RMS_melt_surface_decoy = NaN;
                    RMS_melt_sublithosphere_decoy = NaN;
                    RMS_melt_depths_decoy = NaN;
                else
                    RMS_melt_surface_decoy = [];
                    RMS_melt_sublithosphere_decoy = [];
                    RMS_melt_depths_decoy = [];
                end


                for ttd= 1:ite_depth
                    filename_geofeatures = ['Geofeatures_' model_title '.csv'];
                    fullFilePath_geofeatures = fullfile(path_model_output, filename_geofeatures);
                   
                    % Assuming RMS_melt_surface, RMS_melt_sublithosphere, and RMS_melt_depths 
                    % are defined and ttd is your index or loop variable
                    
                    % Check if any of the values are NaN and replace with 0 if true
 

                    if isnan(RMS_melt_surface_decoy) || isnan(RMS_melt_sublithosphere_decoy) || isnan(RMS_melt_depths_decoy)
                        if isnan(RMS_melt_surface_decoy)
                            RMS_melt_surface = 0;
                        end
                        if isnan(RMS_melt_sublithosphere_decoy)
                            RMS_melt_sublithosphere = 0;
                        end
                        if isnan(RMS_melt_depths_decoy)
                            RMS_melt_depths = 0;
                        end
                    else
                        RMS_melt_surface = RMS_melt_surface(ttd);
                        RMS_melt_sublithosphere = RMS_melt_sublithosphere(ttd);
                        RMS_melt_depths = RMS_melt_depths(ttd);
                    end

                    % Check if the file already exists
                    if exist(fullFilePath_geofeatures, 'file')
                        % File exists, load the existing table
                        geofeatures = readtable(fullFilePath_geofeatures);

                        % Add new values to the existing table
                        newRow = table(time_for_output_structures, numPlumes(ttd), numSubductions(ttd), numSubduction_connected(ttd), current_depth(ttd)./1e3, Maxferet_mean(ttd), Maxferet_max(ttd), Maxferet_min(ttd), nferet(ttd), V_rms_surface, V_rms_continents,RMS_melt_surface,RMS_melt_sublithosphere,RMS_melt_depths, ...
                            'VariableNames', {'Time', 'NumPlumes', 'NumSubductions', 'NumSubduction_connected', 'Depths', 'Maxferet_mean', 'Maxferet_max', 'Maxferet_in', 'nferet', 'V_rms_surface', 'V_rms_continents','surface_rms_melt','sublithosphere_rms_melt','depths_rms_melt'});

                        geofeatures = vertcat(geofeatures, newRow);
                    else
                        % File doesn't exist, create a new table
                        geofeatures = table(time_for_output_structures, numPlumes(ttd), numSubductions(ttd), numSubduction_connected(ttd), current_depth(ttd)./1e3, Maxferet_mean(ttd), Maxferet_max(ttd), Maxferet_min(ttd), nferet(ttd), V_rms_surface, V_rms_continents,RMS_melt_surface,RMS_melt_sublithosphere,RMS_melt_depths,...
                            'VariableNames', {'Time', 'NumPlumes', 'NumSubductions', 'NumSubduction_connected', 'Depths', 'Maxferet_mean', 'Maxferet_max', 'Maxferet_in', 'nferet', 'V_rms_surface', 'V_rms_continents','surface_rms_melt','sublithosphere_rms_melt','depths_rms_melt'});
                    end

                    %         % Display the updated table
                    %     filename_geofeatures = ['Geofeatures_' model_title '.csv'];
                    writetable(geofeatures, fullFilePath_geofeatures);
                end
                disp(geofeatures);
                fprintf('Table saved to %s successfully.\n', filename_geofeatures);
            end

            if any(strcmp(additional_postprocesses, 'oceanic_age_statistics'))

                filename_oceans = ['Geo_oceans_' model_title '.csv'];
                fullFilePath_geooceans = fullfile(path_model_output, filename_oceans);
                oceanics_variable = [time_for_output_structures smoothed_oceanic_age_density]';
                newRow_oceans = table(oceanics_variable, 'VariableNames', {'Time_with_Pdf_oceanic_age'});

                % Check if the file already exists
                if exist(fullFilePath_geooceans, 'file')
                    % File exists, load the existing table
                    geo_oceans = readtable(fullFilePath_geooceans);

                    % Append the new row to the existing table
                    geo_oceans = vertcat(geo_oceans, newRow_oceans);

                    % Write the updated table back to the file
                    writetable(geo_oceans, fullFilePath_geooceans, 'Delimiter', '\t');
                else
                    % File doesn't exist, create a new table with the new row
                    writetable(newRow_oceans, fullFilePath_geooceans, 'Delimiter', '\t');
                end

                fprintf('Table geo_oceans saved to %s successfully.\n', filename_oceans);

            end
            
        end
    catch
        fprintf('For some reason step %s didn''t work and will be skipped. This can occur if a paraview file extracted is missing or corrupted\n', i);
    end
end
end
