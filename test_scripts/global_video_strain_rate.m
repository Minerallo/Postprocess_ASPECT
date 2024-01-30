clc ; close all; clear all;
addpath(genpath('./Individual_scripts/'))

% path_models = {'/Users/ponsm/Desktop/modelblogin/model/globalscale/Potsprocess/sphere3d/P01a_Pangea_1GPa_Mantle_C40MPa_LR/'}
% '/Users/ponsm/Desktop/Model_sphere/Models_HLRN/R01a_Rodinia_1GPa_Mantle_C40MPa_LR/',
path_models = {'/Users/ponsm/Desktop/Model_sphere/Models_HLRN/S01a_no_continents_C40MPA_TP1700_LR/',
            '/Users/ponsm/Desktop/Model_sphere/Models_HLRN/S01h_no_continents_C10MPA_f005_TP1700_LR/'}
% ,
%     '/Users/ponsm/Desktop/Model_sphere/Models_HLRN/S01h_no_continents_C10MPA_f005_TP1700_LR/',
%     '/Users/ponsm/Desktop/Model_sphere/Models_HLRN/S01a_no_continents_C40MPA_TP1700_LR/'}

% Extracting model titles following 'sphere3d/'
model_titles = cell(size(path_models));
for t = 1:numel(path_models)
    splitted = strsplit(path_models{t}, '/');
    idx = find(contains(splitted, 'Models_HLRN'));
    %         'sphere3d'));
    if numel(idx) > 0
        model_titles{t} = splitted{idx+1};
    end
end


for m = 1:numel(path_models)
    path_model = cell2mat(path_models(m));
    make_video='true';
    resample_dynamic_topography = 10;

    %     statistic_parameters = {'Time'};
    %     % Get the data and corresponding statistic numbers and indices
    %     [data_stats, stats_number, header, stat_indices] = get_statistics(path_model, statistic_parameters);
    %
    %     % Statistics
    %     time_vec = data_stats(:, 2)./1e6;
    %     time=num2str(time_vec, '%05.2f Myr');
    %
    %     % find how many zeros for the initial refinement steps
    %     numbers_of_IRS=nnz(~time_vec);

    file_surface = dir(fullfile(path_model, 'surface_*'));
    filenames_surface = {file_surface.name};
    num_files_surface = numel(filenames_surface);

    file_lithosphere = dir(fullfile(path_model, 'lithosphere_*'));
    filenames_lithosphere = {file_lithosphere.name};
    num_files_lithosphere = numel(filenames_lithosphere);

    file_depths = dir(fullfile(path_model, 'depths_*'));
    filenames_depths = {file_depths.name};
    num_files_depths = numel(filenames_depths);

    % Specify the desired frame size for the video
    %     frameSize = [2000, 1000];

    % Create the figure with the specified size
    %     figure('WindowState', 'maximized','Position', [100 100 frameSize]);
    %     scatterSize =70;

    if strcmp(make_video, 'true')
        % Specify the video filename
        videoFilename = [model_titles{m} '.mp4'];

        % Create a VideoWriter object
        videoWriter = VideoWriter(videoFilename, 'MPEG-4');
        videoWriter.FrameRate = 4; % Specify the frame rate (adjust as needed)
        open(videoWriter);
    end

    count_ite = 0;
    Reference_mean_correction_steps=1;
    %     for i = 1:resample_dynamic_topography:num_files_surface

    %     We want the initial position of the continent to be plotted on top of the maps
    surface_init = import_surface_sphere_strain_rate_only(fullfile(path_model, filenames_surface{1}));
    
    % Sort the data in the desired order
    x_surface_init = surface_init.Points0;
    y_surface_init = surface_init.Points1;
    z_surface_init = surface_init.Points2;

    % Calculate the azimuthal angle
    theta_surface_init = atan2(y_surface_init, x_surface_init);
    [theta_sorted_init, theta_index_init] = sort(theta_surface_init);
    x_sorted_surface_init = x_surface_init(theta_index_init);
    y_sorted_surface_init = y_surface_init(theta_index_init);
    z_sorted_surface_init = z_surface_init(theta_index_init);
    surface_init_continent_sorted = surface_init.continent(theta_index_init);
    % Convert Cartesian coordinates to longitude and latitude
    longitude_surface_init = atan2(y_sorted_surface_init, x_sorted_surface_init) * 180 / pi;
    latitude_surface_init = asin(z_sorted_surface_init ./ sqrt(x_sorted_surface_init.^2 + y_sorted_surface_init.^2 + z_sorted_surface_init.^2)) * 180 / pi;
    % Create a grid of longitude and latitude values
    [Xeq, Yeq] = meshgrid(-180:1:180, -90:1:90);
    continents_init = griddata(longitude_surface_init, latitude_surface_init, surface_init_continent_sorted, Xeq, Yeq);
    %         continents_init = continents_init>0.5;
    %         idx=isnan(continents_init);
    %         continents_init(idx)=0;
    %         continents_init = 1 - continents_init;
    %         continents_init(continents_init == 0) = NaN;
    % % %         We binarize the continents for their initial position
    % %         continents_init_binary = imbinarize(continents_init);
    % %
    % %         continent_init_boundaries = bwboundaries(continents_init_binary);
    % %         number_of_init_continent_boundaries = size(continent_init_boundaries, 1);
    % %
    %         continents_init_contour = imcontour(continents_binary);
    %         continents_init_contour(continents_init_contour>360)
%     for i =99
        for i = 1:num_files_surface

        count_ite = count_ite+1;
        fprintf('Processing surface file %d of %d\n', i, num_files_surface)
        file_index = i;

        % Import the data
        surface = import_surface_sphere_strain_rate_only_two(fullfile(path_model, filenames_surface{file_index}));
%         lithosphere = import_lithosphere_sphere(fullfile(path_model, filenames_lithosphere{file_index}));
%         depths = import_depths_sphere(fullfile(path_model, filenames_depths{file_index}));

        %         filename_extract = filenames{file_index};
        %         parts = split(filename_extract, '_');
        %         number_step = str2double(parts{2});
        time_str = sprintf('Time: %.2f Myr', surface.Time(1)./1e6);
        % %05.2f Myr

        % Define the downsampling factor
        downsample_factor = 1;

        % Sort the data in the desired order
        x_surface = surface.Points0;
        y_surface = surface.Points1;
        z_surface = surface.Points2;

%         x_lithosphere = lithosphere.Points0;
%         y_lithosphere = lithosphere.Points1;
%         z_lithosphere = lithosphere.Points2;

%         x_depths = depths.Points0;
%         y_depths = depths.Points1;
%         z_depths = depths.Points2;

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
        surface_topography_sorted = surface_topography(theta_index);
        surface_strain_rate_sorted = surface.strain_rate(theta_index);
        surface_continent_sorted = surface.continent(theta_index);
        surface_v_phi_sorted = surface.v_phi(theta_index);
        surface_v_theta_sorted = surface.v_theta(theta_index);
        surface_v_r_sorted = surface.v_r(theta_index);
%         surface_velocity0_sorted = surface.velocity0(theta_index);
%         surface_velocity1_sorted = surface.velocity1(theta_index);
%         surface_velocity2_sorted = surface.velocity2(theta_index);
% 
%         theta_lithosphere = atan2(y_lithosphere, x_lithosphere);
%         [theta_sorted_lithosphere, theta_index_lithosphere] = sort(theta_lithosphere);
%         x_sorted_lithosphere = x_lithosphere(theta_index_lithosphere);
%         y_sorted_lithosphere = y_lithosphere(theta_index_lithosphere);
%         z_sorted_lithosphere = z_lithosphere(theta_index_lithosphere);
%         divergence_lithosphere_sorted = lithosphere.Divergence(theta_index_lithosphere);
%         oceanic_age_lithosphere_sorted = lithosphere.oceanic_age(theta_index_lithosphere);

%         theta_depths = atan2(y_depths, x_depths);
%         [theta_sorted_depths, theta_index_depths] = sort(theta_depths);
%         x_sorted_depths = x_depths(theta_index_depths);
%         y_sorted_depths = y_depths(theta_index_depths);
%         z_sorted_depths = z_depths(theta_index_depths);
%         depths_sorted_depths = depths.depth(theta_index_depths);
%         unique_depths = unique(depths_sorted_depths);
%         depths_nonadiabatic_temperature_sorted = depths.nonadiabatic_temperature(theta_index_depths);
%         depths_v_phi_sorted = depths.v_phi(theta_index_depths);
%         depths_v_theta_sorted = depths.v_theta(theta_index_depths);
%         depths_v_r_sorted = depths.v_r(theta_index_depths);
%         depths_velocity0_sorted = depths.velocity0(theta_index_depths);
%         depths_velocity1_sorted = depths.velocity1(theta_index_depths);
%         depths_velocity2_sorted = depths.velocity2(theta_index_depths);
% 
%         for u = 1:numel(unique_depths)
%             idx_depths = find(depths_sorted_depths == unique_depths(u));
%             x_sorted_depths_map=x_sorted_depths(idx_depths) ;
%             y_sorted_depths_map=y_sorted_depths(idx_depths);
%             z_sorted_depths_map=z_sorted_depths(idx_depths);
%             x_depths_downsampled{u} = x_sorted_depths_map(1:downsample_factor:end);
%             y_depths_downsampled{u} = y_sorted_depths_map(1:downsample_factor:end);
%             z_depths_downsampled{u} = z_sorted_depths_map(1:downsample_factor:end);
%             non_adiabT_map{u}=depths_nonadiabatic_temperature_sorted(idx_depths);
%             v_theta_map{u}=depths_v_theta_sorted(idx_depths);
%             v_phi_map{u}=depths_v_phi_sorted(idx_depths);
%             v_r_map{u}=depths_v_r_sorted(idx_depths);
%             velocity0_map{u}=depths_velocity0_sorted(idx_depths);
%             velocity1_map{u}=depths_velocity1_sorted(idx_depths);
%             velocity2_map{u}=depths_velocity2_sorted(idx_depths);
%         end

        %         final_surface_topography_sorted = surface_topography_sorted - ((surface_topography_sorted < 2400) .* (1000 / 3300) .* abs(surface_topography_sorted - 2400)) - 2400;
        final_surface_topography_sorted=surface_topography_sorted;


        x_surface_downsampled = x_sorted_surface(1:downsample_factor:end);
        y_surface_downsampled = y_sorted_surface(1:downsample_factor:end);
        z_surface_downsampled = z_sorted_surface(1:downsample_factor:end);
        surface_topography_surface_downsampled = final_surface_topography_sorted(1:downsample_factor:end);
% 
%         x_lithosphere_downsampled = x_sorted_lithosphere(1:downsample_factor:end);
%         y_lithosphere_downsampled = y_sorted_lithosphere(1:downsample_factor:end);
%         z_lithosphere_downsampled = z_sorted_lithosphere(1:downsample_factor:end);

% 
%         min_surface_topography = min(surface_topography_surface_downsampled);
%         max_surface_topography = max(surface_topography_surface_downsampled);
%         mean_surface_topography(count_ite) = mean(surface_topography_surface_downsampled);
% 
%         %        Correction for mean topography subsiding
%         if(count_ite<=Reference_mean_correction_steps)
%             diff_mean_surface_topography(count_ite)=0;
%             disp('test1')
%         elseif (mean_surface_topography(count_ite)<mean_surface_topography(count_ite-1))
%             diff_mean_surface_topography(count_ite)=mean_surface_topography(Reference_mean_correction_steps)-mean_surface_topography(count_ite);
%             disp('test2')
%         else
%             diff_mean_surface_topography(count_ite)=0;
%             disp('test3')
%         end

        % Convert Cartesian coordinates to longitude and latitude
        longitude_surface_downsampled = atan2(y_surface_downsampled, x_surface_downsampled) * 180 / pi;
        latitude_surface_downsampled = asin(z_surface_downsampled ./ sqrt(x_surface_downsampled.^2 + y_surface_downsampled.^2 + z_surface_downsampled.^2)) * 180 / pi;

%         longitude_lithosphere_downsampled = atan2(y_lithosphere_downsampled, x_lithosphere_downsampled) * 180 / pi;
%         latitude_lithosphere_downsampled = asin(z_lithosphere_downsampled ./ sqrt(x_lithosphere_downsampled.^2 + y_lithosphere_downsampled.^2 + z_lithosphere_downsampled.^2)) * 180 / pi;



        % Plot the data on a surface
        %     h2=subplot(1, 2, 2);
        % Create a grid of longitude and latitude values
        [Xeq, Yeq] = meshgrid(-180:1:180, -90:1:90);
        [Xeq_lr, Yeq_lr] = meshgrid(-180:10:180, -90:10:90);
        % Reshape the surface topography data to match the grid size
%         Topo = griddata(longitude_surface_downsampled, latitude_surface_downsampled, surface_topography_surface_downsampled, Xeq, Yeq);
%         Topo = Topo(1:181,1:361)+diff_mean_surface_topography(count_ite);
        V_phi_vec_surface = griddata(longitude_surface_downsampled, latitude_surface_downsampled,surface_v_phi_sorted, Xeq_lr, Yeq_lr);
        V_theta_vec_surface = griddata(longitude_surface_downsampled, latitude_surface_downsampled,surface_v_theta_sorted, Xeq_lr, Yeq_lr);
%         V_phi_surface = griddata(longitude_surface_downsampled, latitude_surface_downsampled,surface_v_phi_sorted, Xeq, Yeq);
%         V_theta_surface = griddata(longitude_surface_downsampled, latitude_surface_downsampled,surface_v_theta_sorted, Xeq, Yeq);
%         V_r_surface = griddata(longitude_surface_downsampled, latitude_surface_downsampled,surface_v_r_sorted, Xeq, Yeq);
%         V_magnitude_surface = sqrt(V_r_surface.^2+V_phi_surface.^2+V_theta_surface.^2);


        strain_rate = griddata(longitude_surface_downsampled, latitude_surface_downsampled, surface_strain_rate_sorted, Xeq, Yeq);
        continents = griddata(longitude_surface_downsampled, latitude_surface_downsampled, surface_continent_sorted, Xeq, Yeq);
%         continents_plot = griddata(longitude_surface_downsampled, latitude_surface_downsampled, surface_continent_sorted, Xeq, Yeq);

%         divergence = griddata(longitude_lithosphere_downsampled, latitude_lithosphere_downsampled, divergence_lithosphere_sorted, Xeq, Yeq);
%         oceanic_age = griddata(longitude_lithosphere_downsampled, latitude_lithosphere_downsampled, oceanic_age_lithosphere_sorted, Xeq, Yeq);
% 
%         for uu = 1:numel(unique_depths)
%             longitude_depths_downsampled = atan2(y_depths_downsampled{uu}, x_depths_downsampled{uu}) * 180 / pi;
%             latitude_depths_downsampled = asin(z_depths_downsampled{uu} ./ sqrt(x_depths_downsampled{uu}.^2 + y_depths_downsampled{uu}.^2 + z_depths_downsampled{uu}.^2)) * 180 / pi;
%             non_adiabT_depths{uu}=griddata(longitude_depths_downsampled, latitude_depths_downsampled, non_adiabT_map{uu}, Xeq, Yeq);
%             V_theta_vec_depths{uu}=griddata(longitude_depths_downsampled, latitude_depths_downsampled, v_theta_map{uu}, Xeq_lr, Yeq_lr);
%             V_phi_vec_depths{uu}=griddata(longitude_depths_downsampled, latitude_depths_downsampled, v_phi_map{uu}, Xeq_lr, Yeq_lr);
%             V_theta_depths{uu}=griddata(longitude_depths_downsampled, latitude_depths_downsampled, v_theta_map{uu}, Xeq, Yeq);
%             V_phi_depths{uu}=griddata(longitude_depths_downsampled, latitude_depths_downsampled, v_phi_map{uu}, Xeq, Yeq);
%             V_r_depths{uu}=griddata(longitude_depths_downsampled, latitude_depths_downsampled, v_r_map{uu}, Xeq, Yeq);
%             Velocity0_depths{uu}=griddata(longitude_depths_downsampled, latitude_depths_downsampled, velocity0_map{uu}, Xeq, Yeq);
%             Velocity1_depths{uu}=griddata(longitude_depths_downsampled, latitude_depths_downsampled, velocity1_map{uu}, Xeq, Yeq);
%             Velocity2_depths{uu}=griddata(longitude_depths_downsampled, latitude_depths_downsampled, velocity2_map{uu}, Xeq, Yeq);
%             V_magnitude_depths{uu}=sqrt(V_r_depths{uu}.^2+V_phi_depths{uu}.^2+V_theta_depths{uu}.^2);
%         end

        %         We do it also for the lithosphere as the coordinate may slightly
        %         change when extracting from paraview which may lead to some
        %         visualisation artefact
        %         continents_lithosphere = griddata(longitude_lithosphere_downsampled, latitude_lithosphere_downsampled, surface_continent_sorted, Xeq, Yeq);


        %       Track continents with a thrshold of 50%
%         continents = continents>0.3;
%         idx=isnan(continents);
%         continents(idx)=0;
%         continents = 1 - continents; %get ocean =1\

        %         continents_lithosphere = continents_lithosphere>0.5;
        %         idx=isnan(continents_lithosphere);
        %         continents_lithosphere(idx)=0;
        %         continents_lithosphere = 1 - continents_lithosphere; %get ocean =1




        %         Bathymetric correction, when there is continent , continents = 0
        %         Topo = Topo - 2400;
        %         Topo = Topo - (continents.*(abs(Topo)* (1000 / 3300)));


        %     mapshow(Xeq,Yeq,Topo,"DisplayType","surface")
        %     surf(Xeq, Yeq, Topo, 'EdgeColor', 'none');
        %     colorbar;
        %     xlabel('Longitude');
        %     ylabel('Latitude');
        %     zlabel('Surface Topography');
        %     title(h2, {'Surface Topography', time_str});c=colorbar;crameri('oleron');caxis([min_surface_topography min_surface_topography+15000]);set(gcf,'color','w'); view(0,90);c.Label.String= "Elevations [m]";

        %         plot topography only
        %         figure(m);
        %         set(gcf, 'Units', 'Normalized', 'OuterPosition', [0, 0.04, 1, 0.96]);
        %         clf;
        %         hold on;
        %         axesm mercator
        %         geoshow(Yeq,Xeq,Topo,'DisplayType', 'surface')
        %         ax = gca;
        %         setm(ax,"MapProjection","robinson")
        %         colorbar;
        %         %     xlabel('Longitude');
        %         %     ylabel('Latitude');
        %         zlabel('Surface Topography');
        %         title({'Surface Topography', time_str});c=colorbar;crameri('oleron');set(gcf,'color','w'); view(0,90);c.Label.String= "Elevations [m]";
        %         caxis([-6000 6000])
        % %         caxis([mean_surface_topography-abs(max_surface_topography-min_surface_topography)/2 mean_surface_topography+abs(max_surface_topography-min_surface_topography)/2])
        %         framem
        %         gridm

        %         figure(m);
        %         set(gcf, 'Units', 'Normalized', 'OuterPosition', [0, 0.04, 1, 0.96]);
        %         clf;
        %         hold on;
        %         axesm mercator
        %         geoshow(Yeq,Xeq,log10(strain_rate),'DisplayType', 'surface')
        %         ax = gca;
        %         setm(ax,"MapProjection","robinson")
        %         colorbar;
        %         %     xlabel('Longitude');
        %         %     ylabel('Latitude');
        %         zlabel('Surface');
        %         title({'Surface', time_str});c=colorbar;crameri('-roma');set(gcf,'color','w'); view(0,90);c.Label.String= "Elevations [m]";
        % %         caxis([-6000 6000])
        % %         caxis([mean_surface_topography-abs(max_surface_topography-min_surface_topography)/2 mean_surface_topography+abs(max_surface_topography-min_surface_topography)/2])
        %         framem
        %         gridm
        %         geoshow(Yeq,Xeq,continents,'facealpha',.3)


%         continents(continents == 0) = NaN;
        %         We binarize the continents for their initial position
%         continents_binary = imbinarize(continents);

        % Assuming 'strain_rate' is the variable you want to process
        
        % Identify negative values in 'strain_rate'
        negative_indices = strain_rate < 0;
        strain_rate_corrected=strain_rate;
        % Replace negative values with NaN
        strain_rate_corrected(negative_indices) = NaN;




        % smoothing_data_coefficient = 0.01;
        % smooth_strain_rate = imgaussfilt(strain_rate,smoothing_data_coefficient);
        scaling_velocity_factor_surface=0.2;
        scaling_velocity_factor_depths=0.01;
        CMap_plume = [1 0 0 ;  1 1 1];
        CMap_subduction = [0 0 1 ;  1 1 1];
        CMap_boundaries = [ 0.0471 0.4392 0.2275;  1 1 1];
        CMap_continents = [ 0.4745 0.3647 0.3020;  1 1 1];
        CMap_continents_green = [ 0 0.4 0.2;  1 1 1];

        % Track continents with a thrshold
        continents = continents>0.3;
        idx=isnan(continents);
        continents(idx)=0;
        continents = 1 - continents;
        %         continents(continents == 0) = NaN;
        continents_position2 = ind2rgb(continents + 1, CMap_continents);
        continents_position_green = ind2rgb(continents + 1, CMap_continents_green);

        Plate_boundaries = strain_rate;
        Plate_boundaries = Plate_boundaries>5e-15;
        Plate_boundaries2 = 1-Plate_boundaries;
%         figure();axesm mercator; geoshow(Xeq,Yeq,Plate_boundaries)

        %% Track plumes and subductions over depths

        plumes_depths_tracking=[1,5];
        %select the depth files 1 is the shallowest and get deeper as the value increase
        subduction_depths_tracking=[1,5];
%         displayed_depth = u;
        plumes_non_adiabtic_tracking_temperature = [30,50];%5
        subduction_non_adiabtic_tracking_temperature = [-200,-400];%5



        %% Plot Strain rate
        figure(1);
        set(gcf, 'Units', 'Normalized', 'OuterPosition', [0, 0.04, 1, 0.96]);
        clf;
        set(gcf, 'color', 'w');
        view(0, 90);
        axesm mercator

        geoshow(Yeq, Xeq, log10(strain_rate_corrected), 'DisplayType', 'surface');
        ax1 = gca;
        setm(ax1, 'MapProjection', 'robinson')
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
%         figure();        axesm mercator

        ax2 = gca;
        setm(ax2, 'MapProjection', 'robinson')
        geoshow(Yeq, Xeq, continents_position_green,'facealpha', 0.3);
        colormap(ax2);freezeColors;

        contourm(Yeq,Xeq, continents_init, 0.3, 'w','LineWidth',2);
        quivermc(Yeq_lr, Xeq_lr, V_theta_vec_surface, V_phi_vec_surface,'color','k','linewidth',1,'units','Velocity (m.yr^-^1)','reference',scaling_velocity_factor_surface);
        ax2.SortMethod='childorder';


      

        % Tile all open figures
%         tilefigs;

    

        if strcmp(make_video, 'true')
            % Capture the current frame
            frame = getframe(gcf);

            % Write the frame to the video file
            writeVideo(videoWriter, frame);
        end

    end

    if strcmp(make_video, 'true')
        % Close the video file
        close(videoWriter);
    end
end