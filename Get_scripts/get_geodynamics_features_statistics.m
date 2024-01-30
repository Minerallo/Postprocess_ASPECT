function [Geodynamics_structures_distribution] = get_geodynamics_features_statistics(path_model, additional_postprocesses)

allowed_parameters = {'subduction_and_plume_statistics','oceanic_age_statistics','continents_RMS'};
% Check if averaged_parameter is allowed
if ~all(ismember(additional_postprocesses, allowed_parameters))
    error('Invalid parameter selected. The allowed parameters are: %s', strjoin(allowed_parameters, ', '));
end

Geodynamics_structures_distribution = table();

% Split the path into parts
parts = strsplit(path_model, '/');

% Extract the last part
model_title = parts{end-1};

%initiate these variable for further checking
num_files_surface=[];
num_files_lithosphere=[];
num_files_depths=[];
files_required = [];
files_required_2 = [];
files_required_3 = [];
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
%Required files for each scheme
if any(strcmp(additional_postprocesses, 'subduction_and_plume_statistics'))
    files_required_1 = {'surface','depths'};
end
if any(strcmp(additional_postprocesses, 'oceanic_age_statistics'))
    files_required_2 = {'lithosphere'};
end
if any(strcmp(additional_postprocesses, 'continents_RMS'))
    files_required_3 = {'surface'};
end
total_files_required = [files_required_1 files_required_2 files_required_3];
files_requirement = unique(total_files_required);

% Check if subduction_and_plume_number is true
if any(strcmp(files_requirement, 'surface'))
    file_surface = dir(fullfile(path_model, 'surface_*'));
    filenames_surface = {file_surface.name};
    num_files_surface = numel(filenames_surface);
    disp(['Found ' num2str(num_files_surface) ' surface files.']);

end

if any(strcmp(files_requirement, 'depths'))
    file_depths = dir(fullfile(path_model, 'depths_*'));
    filenames_depths = {file_depths.name};
    num_files_depths = numel(filenames_depths);
    disp(['Found ' num2str(num_files_depths) ' depths files.']);

end


if any(strcmp(files_requirement, 'lithosphere'))
    file_lithosphere = dir(fullfile(path_model, 'lithosphere_*'));
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

% Check if the number of files for surface and depths is different
if num_files_surface ~= num_files_depths
    fprintf('Number of surface files (%d) is different from the number of depth files (%d), the satistics will only be counted until file number (%d).', num_files_surface, num_files_depths,min_num_files);
end

count_ite = 0;
Reference_mean_correction_steps=1;

% Create a grid of longitude and latitude values
[Xeq, Yeq] = meshgrid(-180:1:180, -90:1:90);

for i = 1:min_num_files

    count_ite = count_ite+1;
    fprintf('Processing additional statistics files %d of %d\n', i, min_num_files)
    file_index = i;

    % Import only the ncessary data
    if  any(strcmp(files_requirement, 'surface'))
        surface = import_surface_sphere(fullfile(path_model, filenames_surface{file_index}));
    end
    if  any(strcmp(files_requirement, 'depths'))
        depths = import_depths_sphere(fullfile(path_model, filenames_depths{file_index}));
    end
    if any(strcmp(files_requirement, 'lithosphere'))
        lithosphere = import_lithosphere_sphere(fullfile(path_model, filenames_lithosphere{file_index}));
    end

    % Lets define where the time should be read
    if any(strcmp(files_requirement,'surface'))
        time_str = sprintf('Time: %.2f Myr', surface.Time(1) / 1e6);
        time_for_output_structures = surface.Time(1)./1e6;
    elseif any(strcmp(files_requirement,'lithosphere'))
        time_str = sprintf('Time: %.2f Myr', lithosphere.Time(1) / 1e6);
        time_for_output_structures = lithosphere.Time(1)./1e6;
    elseif any(strcmp(files_requirement,'depths'))
        time_str = sprintf('Time: %.2f Myr', depths.Time(1) / 1e6);
        time_for_output_structures = depths.Time(1)./1e6;
    end
    % Define the downsampling factor
    downsample_factor = 1;

    if  any(strcmp(files_requirement, 'surface'))

        % Sort the data in the desired order
        x_surface = surface.Points0;
        y_surface = surface.Points1;
        z_surface = surface.Points2;

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
        surface_velocity0_sorted = surface.velocity0(theta_index);
        surface_velocity1_sorted = surface.velocity1(theta_index);
        surface_velocity2_sorted = surface.velocity2(theta_index);
        surface_topography_sorted = surface_topography(theta_index);
        surface_strain_rate_sorted = surface.strain_rate(theta_index);
        surface_continent_sorted = surface.continent(theta_index);
        final_surface_topography_sorted=surface_topography_sorted;


        x_surface_downsampled = x_sorted_surface(1:downsample_factor:end);
        y_surface_downsampled = y_sorted_surface(1:downsample_factor:end);
        z_surface_downsampled = z_sorted_surface(1:downsample_factor:end);
        surface_topography_surface_downsampled = final_surface_topography_sorted(1:downsample_factor:end);
        % Convert Cartesian coordinates to longitude and latitude
        longitude_surface_downsampled = atan2(y_surface_downsampled, x_surface_downsampled) * 180 / pi;
        latitude_surface_downsampled = asin(z_surface_downsampled ./ sqrt(x_surface_downsampled.^2 + y_surface_downsampled.^2 + z_surface_downsampled.^2)) * 180 / pi;
    end

    if any(strcmp(files_requirement, 'depths'))
        x_depths = depths.Points0;
        y_depths = depths.Points1;
        z_depths = depths.Points2;

        theta_depths = atan2(y_depths, x_depths);
        [theta_sorted_depths, theta_index_depths] = sort(theta_depths);
        x_sorted_depths = x_depths(theta_index_depths);
        y_sorted_depths = y_depths(theta_index_depths);
        z_sorted_depths = z_depths(theta_index_depths);
        depths_sorted_depths = depths.depth(theta_index_depths);
        unique_depths = unique(depths_sorted_depths);
        depths_nonadiabatic_temperature_sorted = depths.nonadiabatic_temperature(theta_index_depths);


        for u = 1:numel(unique_depths)
            idx_depths = find(depths_sorted_depths == unique_depths(u));
            x_sorted_depths_map=x_sorted_depths(idx_depths) ;
            y_sorted_depths_map=y_sorted_depths(idx_depths);
            z_sorted_depths_map=z_sorted_depths(idx_depths);
            non_adiabT_map{u}=depths_nonadiabatic_temperature_sorted(idx_depths);
            x_depths_downsampled{u} = x_sorted_depths_map(1:downsample_factor:end);
            y_depths_downsampled{u} = y_sorted_depths_map(1:downsample_factor:end);
            z_depths_downsampled{u} = z_sorted_depths_map(1:downsample_factor:end);
        end
    end

    %Lets set griddata separately so we only call it only if needed
    if  any(strcmp(additional_postprocesses, 'subduction_and_plume_statistics'))

        min_surface_topography = min(surface_topography_surface_downsampled);
        max_surface_topography = max(surface_topography_surface_downsampled);
        mean_surface_topography(count_ite) = mean(surface_topography_surface_downsampled);

        %        Correction for mean topography subsiding
        if(count_ite<=Reference_mean_correction_steps)
            diff_mean_surface_topography(count_ite)=0;
            %             disp('test1')
        elseif (mean_surface_topography(count_ite)<mean_surface_topography(count_ite-1))
            diff_mean_surface_topography(count_ite)=mean_surface_topography(Reference_mean_correction_steps)-mean_surface_topography(count_ite);
            %             disp('test2')
        else
            diff_mean_surface_topography(count_ite)=0;
            %             disp('test3')
        end

        % Reshape the surface topography data to match the grid size
        Topo = griddata(longitude_surface_downsampled, latitude_surface_downsampled, surface_topography_surface_downsampled, Xeq, Yeq);
        Topo = Topo(1:181,1:361)+diff_mean_surface_topography(count_ite);

        for uu = 1:numel(unique_depths)
            longitude_depths_downsampled = atan2(y_depths_downsampled{uu}, x_depths_downsampled{uu}) * 180 / pi;
            latitude_depths_downsampled = asin(z_depths_downsampled{uu} ./ sqrt(x_depths_downsampled{uu}.^2 + y_depths_downsampled{uu}.^2 + z_depths_downsampled{uu}.^2)) * 180 / pi;
            non_adiabT_depths{uu}=griddata(longitude_depths_downsampled, latitude_depths_downsampled, non_adiabT_map{uu}, Xeq, Yeq);
        end
    end

    if any(strcmp(files_requirement, 'lithosphere'))
        x_lithosphere = lithosphere.Points0;
        y_lithosphere = lithosphere.Points1;
        z_lithosphere = lithosphere.Points2;

        theta_lithosphere = atan2(y_lithosphere, x_lithosphere);
        [theta_sorted_lithosphere, theta_index_lithosphere] = sort(theta_lithosphere);
        x_sorted_lithosphere = x_lithosphere(theta_index_lithosphere);
        y_sorted_lithosphere = y_lithosphere(theta_index_lithosphere);
        z_sorted_lithosphere = z_lithosphere(theta_index_lithosphere);
        x_sorted_lithosphere = x_lithosphere(theta_index_lithosphere);
        y_sorted_lithosphere = y_lithosphere(theta_index_lithosphere);
        z_sorted_lithosphere = z_lithosphere(theta_index_lithosphere);
        oceanic_age_lithosphere_sorted = lithosphere.oceanic_age(theta_index_lithosphere);
        lithosphere_continent_sorted = lithosphere.continent(theta_index_lithosphere);

        x_lithosphere_downsampled = x_sorted_lithosphere(1:downsample_factor:end);
        y_lithosphere_downsampled = y_sorted_lithosphere(1:downsample_factor:end);
        z_lithosphere_downsampled = z_sorted_lithosphere(1:downsample_factor:end);

        longitude_lithosphere_downsampled = atan2(y_lithosphere_downsampled, x_lithosphere_downsampled) * 180 / pi;
        latitude_lithosphere_downsampled = asin(z_lithosphere_downsampled ./ sqrt(x_lithosphere_downsampled.^2 + y_lithosphere_downsampled.^2 + z_lithosphere_downsampled.^2)) * 180 / pi;
    end

    %Lets set griddata separately so we only call it only if needed
    if any(strcmp(additional_postprocesses, 'oceanic_age_statistics'))
        oceanic_age = griddata(longitude_lithosphere_downsampled, latitude_lithosphere_downsampled, oceanic_age_lithosphere_sorted, Xeq, Yeq);
    end

    if any(strcmp(files_requirement, 'surface'))
        continents = griddata(longitude_surface_downsampled, latitude_surface_downsampled, surface_continent_sorted, Xeq, Yeq);
    elseif  ~any(strcmp(files_requirement, 'surface'))
        if any(strcmp(files_requirement, 'lithosphere'))
            continents = griddata(longitude_lithosphere_downsampled, latitude_lithosphere_downsampled, lithosphere_continent_sorted, Xeq, Yeq);
        end
    end


    %% Statistics place
    if any(strcmp(additional_postprocesses, 'subduction_and_plume_statistics'))
        %% Plot topography corrected

        %continents positions
        continents_threshold = continents>0.3;
        idx_is_nan_to_zero=isnan(continents_threshold);
        continents_threshold(idx_is_nan_to_zero)=0;
        oceans = 1 - continents_threshold;
        %with resolution of Xeq and Yeq of 1deg from griddata
        idx_oceans=find(oceans == 1);
        idx_cont=find(oceans == 0);

        % Correct the topography for the water load in oceans
        smooth_Topo = imgaussfilt(Topo,0.3);
        Topo_corrected=smooth_Topo ;
        Topo_corrected(idx_oceans) = Topo_corrected(idx_oceans) - (Topo_corrected(idx_oceans)<0).* (1000 / 3300).* abs(Topo_corrected(idx_oceans));

        %% Track plumes and subductions over depths

        plumes_depths_tracking=[1,5];
        %select the depth files 1 is the shallowest and get deeper as the value increase
        subduction_depths_tracking=[1,5];
        displayed_depth = u;
        plumes_non_adiabtic_tracking_temperature = [30,150];%5
        subduction_non_adiabtic_tracking_temperature = [-200,-400];%5

        for tt= 1:numel(plumes_depths_tracking)
            %%

            current_depth=unique_depths(plumes_depths_tracking(tt));
            plumes_position = non_adiabT_depths{plumes_depths_tracking(tt)}>plumes_non_adiabtic_tracking_temperature(tt);
            subductions_position = non_adiabT_depths{subduction_depths_tracking(tt)}<subduction_non_adiabtic_tracking_temperature(tt);
            plumes_positions2=1-plumes_position;
            subductions_positions2=1-subductions_position;
            %             subductions_position3 = ind2rgb(subductions_positions2 + 1, CMap_subduction);
            %             boundaries_position3 = ind2rgb(Plate_boundaries2+1, CMap_boundaries);
            cmatrix_plume = contourm(Yeq,Xeq, non_adiabT_depths{plumes_depths_tracking(tt)}, [plumes_non_adiabtic_tracking_temperature(tt) plumes_non_adiabtic_tracking_temperature(tt)], 'r','LineWidth',2);
            cmatrix_subduction = contourm(Yeq,Xeq, non_adiabT_depths{subduction_depths_tracking(tt)}, [-90000 subduction_non_adiabtic_tracking_temperature(tt)], 'b','LineWidth',2);
            %            cmatrix_boundaries = contourm(Yeq,Xeq, log10(strain_rate), [5e-15 1], 'm','LineWidth',2);
            %%
            if isempty(cmatrix_plume)
                numPlumes = 0;
            else
                [x_plume, y_plume] = C2xyz(cmatrix_plume);
                polyin_plumes = polyshape(x_plume, y_plume);
                numPlumes = polyin_plumes.NumRegions;
            end
            if isempty(cmatrix_subduction)
                numSubductions = 0;
                numSubduction_connected = 0;
                Maxferet_mean = 0;
                Maxferet_max = 0;
                Maxferet_min = 0;
                nferet = 0;
            else
                [x_subduction, y_subduction] = C2xyz(cmatrix_subduction);
                polyin_subductions = polyshape(x_subduction, y_subduction);
                numSubductions = polyin_subductions.NumRegions;
                It2 = bwmorph(subductions_position,'thin','inf');
                %                 If some cleaning is needed
                %                 It2 = bwmorph(It,'clean');
                Subduction_connectivity= bwconncomp(It2);
                numSubduction_connected = Subduction_connectivity.NumObjects;
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
                Maxferet_mean = mean(distance_feret_km);
                Maxferet_max = max(distance_feret_km);
                Maxferet_min = min(distance_feret_km);
                nferet = numel(distance_feret_km);

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
        end
    end

    if any(strcmp(additional_postprocesses, 'oceanic_age_statistics'))
        %% Plot oceanic age
        %         idx_oceanic_domain = find(oceanic == 0);
        %         oceanic_age_domain = oceanic_age;
        %         oceanic_age_domain(continents_threshold) = NaN;
        % Get the active plate boundaries using a threshold on strain rate

        % Plot similar to Coltice et al 2012, distribution of oceanic age
        %
        % bin_edges_oce_age = 0:1:250;
        %
        % % Create the histogram
        % figure();
        % bin_counts = histcounts(oceanic_age_domain, bin_edges_oce_age);
        % area_per_unit_of_age = bin_counts.*111*111;
        % histogram(oceanic_age_domain, bin_edges_oce_age, 'FaceColor', 'blue', 'EdgeColor', 'black');
        % xlabel('Oceanic Age (My)');
        % ylabel('Frequency');
        % title('Histogram of Oceanic Age');
        % figure();
        % plot(bin_edges_oce_age(2:end),area_per_unit_of_age);xlabel('Oceanic Age (My)');ylabel('Area per unit of age [km^2 yr^-^1]');
        % title('Oceanic age distribution');
    end

    if any(strcmp(additional_postprocesses, 'continents_RMS'))
        %Lets avoid to use griddata, lets calculate the index
        idx_continents_sorted = find(surface_continent_sorted>0.3);
        norm_square_vel_surface = sqrt(surface_velocity0_sorted.^2+surface_velocity1_sorted.^2+surface_velocity2_sorted.^2);
        norm_square_vel_continents = sqrt(surface_velocity0_sorted(idx_continents_sorted).^2+surface_velocity1_sorted(idx_continents_sorted).^2+surface_velocity2_sorted(idx_continents_sorted).^2);
        V_rms_surface = rms(norm_square_vel_surface);
        V_rms_continents = rms(norm_square_vel_continents);
    end

    %% Write everything down
    % We already check in postprocess.m whether a file already exists, so that
    % we plot the data directly from the file and do not rewrite the file.
    % However, if someone wants to add data to the same file again
    % and again plotting can be switched off with postprocess.m
    % and the file could be opened and written again, thanks to the first
    % condition ~exist('Geodynamics_structures_distribution', 'var')

    % Check if the table exists, if not, create a new table
    if ~exist('Geodynamics_structures_distribution', 'var')
        Geodynamics_structures_distribution = table(time_for_output_structures, numPlumes, numSubductions,numSubduction_connected,current_depth./1e3,Maxferet_mean,Maxferet_max,Maxferet_min,nferet,V_rms_surface,V_rms_continents, ...
            'VariableNames', {'Time', 'NumPlumes', 'NumSubductions','NumSubduction_connected','Depths','Maxferet_mean','Maxferet_max','Maxferet_in','nferet','V_rms_surface','V_rms_continents'});
    else
        % Add new values to the existing table
        newRow = table(time_for_output_structures, numPlumes, numSubductions,numSubduction_connected,current_depth./1e3,Maxferet_mean,Maxferet_max,Maxferet_min,nferet,V_rms_surface,V_rms_continents, ...
            'VariableNames', {'Time', 'NumPlumes', 'NumSubductions','NumSubduction_connected','Depths','Maxferet_mean','Maxferet_max','Maxferet_in','nferet','V_rms_surface','V_rms_continents'});

        Geodynamics_structures_distribution = vertcat(Geodynamics_structures_distribution, newRow);
    end

    %         % Display the updated table
    disp(Geodynamics_structures_distribution);

    % Save the updated table to a file (adjust the filename and format as needed)
    writetable(Geodynamics_structures_distribution, ['Geofeatures_' model_title '.csv']);
    %                 ['Distribution_' model_titles{m} '.csv']);

end

end
