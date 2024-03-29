%This file project the data extracted in Paraview for spherical 3D models
%using the Global extract_series.py. 
% The data are projected using a robinson projection
% The file is till being built so sorry for the mess

clc ; close all; clear all;
addpath(genpath('./Individual_scripts/'))

path_models = {'/Users/ponsm/Desktop/Model_sphere/Models_HLRN/P01a_Pangea_1GPa_Mantle_C40MPa_LR/'}

% Extracting model titles
model_titles = cell(size(path_models));
for t = 1:numel(path_models)
% Split the path into parts
parts = strsplit(path_models{t}, '/');
% Extract the last part
model_titles{t} = parts{end-1};
end

Data_to_visualize = {'Topography','Geofeatures','Strain_rate', 'Velocity_magnitude','non_adiabatic_temperature'}; 
make_video='false';

if strcmp(make_video,'false')
end

for m = 1:numel(path_models)
    path_model = cell2mat(path_models(m));
%     resample_dynamic_topography = 10;

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

%     count_ite = 0;
    Reference_mean_correction_steps=2;

    %     We want the initial position of the continent to be plotted on top of the maps
    surface_init = import_surface_sphere(fullfile(path_model, filenames_surface{Reference_mean_correction_steps}));
    % Sort the data in the desired order
    x_surface_init = surface_init.Points0;
    y_surface_init = surface_init.Points1;
    z_surface_init = surface_init.Points2;

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
    surface_init_continent_sorted = surface_init.continent(theta_index_init);
    % Convert Cartesian coordinates to longitude and latitude
    longitude_surface_init = atan2(y_sorted_surface_init, x_sorted_surface_init) * 180 / pi;
    latitude_surface_init = asin(z_sorted_surface_init ./ sqrt(x_sorted_surface_init.^2 + y_sorted_surface_init.^2 + z_sorted_surface_init.^2)) * 180 / pi;
    % Create a grid of longitude and latitude values
    [Xeq, Yeq] = meshgrid(-180:1:180, -90:1:90);
    [Xeq_lr, Yeq_lr] = meshgrid(-180:5:180, -90:5:90);
    continents_init = griddata(longitude_surface_init, latitude_surface_init, surface_init_continent_sorted, Xeq, Yeq);

    mean_surface_init_topography = mean(surface_init_topography_sorted);


    for i =99
        %     for i = 1:num_files_surface

%         count_ite = i-1;
        count_ite = i;
        fprintf('Processing surface file %d of %d\n', i, num_files_surface)
        file_index = i;

        % Import the data
        surface = import_surface_sphere(fullfile(path_model, filenames_surface{file_index}));
        lithosphere = import_lithosphere_sphere(fullfile(path_model, filenames_lithosphere{file_index}));
        depths = import_depths_sphere(fullfile(path_model, filenames_depths{file_index}));

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

        x_lithosphere = lithosphere.Points0;
        y_lithosphere = lithosphere.Points1;
        z_lithosphere = lithosphere.Points2;

        x_depths = depths.Points0;
        y_depths = depths.Points1;
        z_depths = depths.Points2;

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
        surface_velocity0_sorted = surface.velocity0(theta_index);
        surface_velocity1_sorted = surface.velocity1(theta_index);
        surface_velocity2_sorted = surface.velocity2(theta_index);

        theta_lithosphere = atan2(y_lithosphere, x_lithosphere);
        [theta_sorted_lithosphere, theta_index_lithosphere] = sort(theta_lithosphere);
        x_sorted_lithosphere = x_lithosphere(theta_index_lithosphere);
        y_sorted_lithosphere = y_lithosphere(theta_index_lithosphere);
        z_sorted_lithosphere = z_lithosphere(theta_index_lithosphere);
        divergence_lithosphere_sorted = lithosphere.Divergence(theta_index_lithosphere);
        oceanic_age_lithosphere_sorted = lithosphere.oceanic_age(theta_index_lithosphere);

        theta_depths = atan2(y_depths, x_depths);
        [theta_sorted_depths, theta_index_depths] = sort(theta_depths);
        x_sorted_depths = x_depths(theta_index_depths);
        y_sorted_depths = y_depths(theta_index_depths);
        z_sorted_depths = z_depths(theta_index_depths);
        depths_sorted_depths = depths.depth(theta_index_depths);
        unique_depths = unique(depths_sorted_depths);
        depths_nonadiabatic_temperature_sorted = depths.nonadiabatic_temperature(theta_index_depths);
        depths_v_phi_sorted = depths.v_phi(theta_index_depths);
        depths_v_theta_sorted = depths.v_theta(theta_index_depths);
        depths_v_r_sorted = depths.v_r(theta_index_depths);
        depths_velocity0_sorted = depths.velocity0(theta_index_depths);
        depths_velocity1_sorted = depths.velocity1(theta_index_depths);
        depths_velocity2_sorted = depths.velocity2(theta_index_depths);

        for u = 1:numel(unique_depths)
            idx_depths = find(depths_sorted_depths == unique_depths(u));
            x_sorted_depths_map=x_sorted_depths(idx_depths) ;
            y_sorted_depths_map=y_sorted_depths(idx_depths);
            z_sorted_depths_map=z_sorted_depths(idx_depths);
            x_depths_downsampled{u} = x_sorted_depths_map(1:downsample_factor:end);
            y_depths_downsampled{u} = y_sorted_depths_map(1:downsample_factor:end);
            z_depths_downsampled{u} = z_sorted_depths_map(1:downsample_factor:end);
            non_adiabT_map{u}=depths_nonadiabatic_temperature_sorted(idx_depths);
            v_theta_map{u}=depths_v_theta_sorted(idx_depths);
            v_phi_map{u}=depths_v_phi_sorted(idx_depths);
            v_r_map{u}=depths_v_r_sorted(idx_depths);
            velocity0_map{u}=depths_velocity0_sorted(idx_depths);
            velocity1_map{u}=depths_velocity1_sorted(idx_depths);
            velocity2_map{u}=depths_velocity2_sorted(idx_depths);
        end

        %         final_surface_topography_sorted = surface_topography_sorted - ((surface_topography_sorted < 2400) .* (1000 / 3300) .* abs(surface_topography_sorted - 2400)) - 2400;
        final_surface_topography_sorted=surface_topography_sorted;


        x_surface_downsampled = x_sorted_surface(1:downsample_factor:end);
        y_surface_downsampled = y_sorted_surface(1:downsample_factor:end);
        z_surface_downsampled = z_sorted_surface(1:downsample_factor:end);
        surface_topography_surface_downsampled = final_surface_topography_sorted(1:downsample_factor:end);

        x_lithosphere_downsampled = x_sorted_lithosphere(1:downsample_factor:end);
        y_lithosphere_downsampled = y_sorted_lithosphere(1:downsample_factor:end);
        z_lithosphere_downsampled = z_sorted_lithosphere(1:downsample_factor:end);


        min_surface_topography = min(surface_topography_surface_downsampled);
        max_surface_topography = max(surface_topography_surface_downsampled);
        mean_surface_topography(count_ite) = mean(surface_topography_surface_downsampled);

        %        Correction for mean topography subsiding
        if(count_ite<=Reference_mean_correction_steps)
            diff_mean_surface_topography(count_ite)=0;
            disp('test1')
        elseif (mean_surface_topography(count_ite)<mean_surface_topography(count_ite-1))
            diff_mean_surface_topography(count_ite)=mean_surface_init_topography-mean_surface_topography(count_ite);
            disp('test2')
        else
            diff_mean_surface_topography(count_ite)=0;
            disp('test3')
        end

        % Convert Cartesian coordinates to longitude and latitude
        longitude_surface_downsampled = atan2(y_surface_downsampled, x_surface_downsampled) * 180 / pi;
        latitude_surface_downsampled = asin(z_surface_downsampled ./ sqrt(x_surface_downsampled.^2 + y_surface_downsampled.^2 + z_surface_downsampled.^2)) * 180 / pi;

        longitude_lithosphere_downsampled = atan2(y_lithosphere_downsampled, x_lithosphere_downsampled) * 180 / pi;
        latitude_lithosphere_downsampled = asin(z_lithosphere_downsampled ./ sqrt(x_lithosphere_downsampled.^2 + y_lithosphere_downsampled.^2 + z_lithosphere_downsampled.^2)) * 180 / pi;


        % Plot the data on a surface
        %     h2=subplot(1, 2, 2);
        % Reshape the surface topography data to match the grid size
        Topo = griddata(longitude_surface_downsampled, latitude_surface_downsampled, surface_topography_surface_downsampled, Xeq, Yeq);
        Topo = Topo(1:181,1:361)+diff_mean_surface_topography(count_ite);
        V_phi_vec_surface = griddata(longitude_surface_downsampled, latitude_surface_downsampled,surface_v_phi_sorted, Xeq_lr, Yeq_lr);
        V_theta_vec_surface = griddata(longitude_surface_downsampled, latitude_surface_downsampled,surface_v_theta_sorted, Xeq_lr, Yeq_lr);
        V_phi_surface = griddata(longitude_surface_downsampled, latitude_surface_downsampled,surface_v_phi_sorted, Xeq, Yeq);
        V_theta_surface = griddata(longitude_surface_downsampled, latitude_surface_downsampled,surface_v_theta_sorted, Xeq, Yeq);
        V_r_surface = griddata(longitude_surface_downsampled, latitude_surface_downsampled,surface_v_r_sorted, Xeq, Yeq);
        V_magnitude_surface = sqrt(V_r_surface.^2+V_phi_surface.^2+V_theta_surface.^2);

        strain_rate = griddata(longitude_surface_downsampled, latitude_surface_downsampled, surface_strain_rate_sorted, Xeq, Yeq);
        continents = griddata(longitude_surface_downsampled, latitude_surface_downsampled, surface_continent_sorted, Xeq, Yeq);
        continents_plot = griddata(longitude_surface_downsampled, latitude_surface_downsampled, surface_continent_sorted, Xeq, Yeq);

        divergence = griddata(longitude_lithosphere_downsampled, latitude_lithosphere_downsampled, divergence_lithosphere_sorted, Xeq, Yeq);
        oceanic_age = griddata(longitude_lithosphere_downsampled, latitude_lithosphere_downsampled, oceanic_age_lithosphere_sorted, Xeq, Yeq);

        for uu = 1:numel(unique_depths)
            longitude_depths_downsampled = atan2(y_depths_downsampled{uu}, x_depths_downsampled{uu}) * 180 / pi;
            latitude_depths_downsampled = asin(z_depths_downsampled{uu} ./ sqrt(x_depths_downsampled{uu}.^2 + y_depths_downsampled{uu}.^2 + z_depths_downsampled{uu}.^2)) * 180 / pi;
            non_adiabT_depths{uu}=griddata(longitude_depths_downsampled, latitude_depths_downsampled, non_adiabT_map{uu}, Xeq, Yeq);
            V_theta_vec_depths{uu}=griddata(longitude_depths_downsampled, latitude_depths_downsampled, v_theta_map{uu}, Xeq_lr, Yeq_lr);
            V_phi_vec_depths{uu}=griddata(longitude_depths_downsampled, latitude_depths_downsampled, v_phi_map{uu}, Xeq_lr, Yeq_lr);
            V_theta_depths{uu}=griddata(longitude_depths_downsampled, latitude_depths_downsampled, v_theta_map{uu}, Xeq, Yeq);
            V_phi_depths{uu}=griddata(longitude_depths_downsampled, latitude_depths_downsampled, v_phi_map{uu}, Xeq, Yeq);
            V_r_depths{uu}=griddata(longitude_depths_downsampled, latitude_depths_downsampled, v_r_map{uu}, Xeq, Yeq);
            Velocity0_depths{uu}=griddata(longitude_depths_downsampled, latitude_depths_downsampled, velocity0_map{uu}, Xeq, Yeq);
            Velocity1_depths{uu}=griddata(longitude_depths_downsampled, latitude_depths_downsampled, velocity1_map{uu}, Xeq, Yeq);
            Velocity2_depths{uu}=griddata(longitude_depths_downsampled, latitude_depths_downsampled, velocity2_map{uu}, Xeq, Yeq);
            V_magnitude_depths{uu}=sqrt(V_r_depths{uu}.^2+V_phi_depths{uu}.^2+V_theta_depths{uu}.^2);
        end


        % smoothing_data_coefficient = 0.01;
        % smooth_strain_rate = imgaussfilt(strain_rate,smoothing_data_coefficient);
        scaling_velocity_factor_surface=0.2;
        scaling_velocity_factor_depths=0.01;
        CMap_plume = [1 0 0 ;  1 1 1];
        CMap_subduction = [0 0 1 ;  1 1 1];
        CMap_boundaries = [ 0.0471 0.4392 0.2275;  1 1 1];
        CMap_continents = [ 0.4745 0.3647 0.3020;  1 1 1];
        CMap_continents_dark = [ 0 0 0;  1 1 1];

        % Track continents with a thrshold
        continents_threshold = continents>0.3;
        idx_is_nan_to_zero=isnan(continents_threshold);
        continents_threshold(idx_is_nan_to_zero)=0;
        oceans = 1 - continents_threshold;
        %         continents(continents == 0) = NaN;
        continents_reversed_position_for_plot = ind2rgb(oceans + 1, CMap_continents);
        continents_reversed_position_for_plot_dark = ind2rgb(oceans + 1, CMap_continents_dark);

        Plate_boundaries = strain_rate;
        Plate_boundaries = Plate_boundaries>5e-15;
        Plate_boundaries2 = 1-Plate_boundaries;
%         figure();axesm mercator; geoshow(Xeq,Yeq,Plate_boundaries)

%% Plot topography corrected 
 smooth_Topo = imgaussfilt(Topo,0.3);

 Topo_corrected=smooth_Topo ; 
 idx_oceans=find(oceans == 1);
  idx_cont=find(oceans == 0);
Topo_corrected(idx_oceans) = Topo_corrected(idx_oceans) - (Topo_corrected(idx_oceans)<0).* (1000 / 3300).* abs(Topo_corrected(idx_oceans));
% -2400 was removed
%%
Topography_colors = [...
      -11000          10           0         121;
      -10500          26           0         137;
      -10000          38           0         152;
       -9500          27           3         166;
       -9000          16           6         180;
       -8500           5           9         193;
       -8000           0          14         203;
       -7500           0          22         210;
       -7000           0          30         216;
       -6500           0          39         223;
       -6000          12          68         231;
       -5500          26         102         240;
       -5000          19         117         244;
       -4500          14         133         249;
       -4000          21         158         252;
       -3500          30         178         255;
       -3000          43         186         255;
       -2500          55         193         255;
       -2000          65         200         255;
       -1500          79         210         255;
       -1000          94         223         255;
        -500         138         227         255;
          0          16         123          48;
         500         232         214         125;
        1200         163          68           0;
        1800         130          30          30;
        2800         161         161         161;
        3000         206         206         206;
        3500         200         200         200;
        4000         180         180         180;
        4500         160         160         160;
        5000          90          90          90];
        figure(1);
        set(gcf, 'Units', 'Normalized', 'OuterPosition', [0, 0.04, 1, 0.96],'color','w');clf;
        view(0, 90);
        axesm mercator
        framem;
        gridm;
        ax0 = gca;
        setm(ax0, 'MapProjection', 'robinson');
%         geoshow(Yeq, Xeq,Topo_corrected);
        geoimg = geoshow(Yeq, Xeq,Topo_corrected,'DisplayType', 'texturemap'); % save object
        shadem(-14.5,[210 75]) 
        geoimg.AlphaDataMapping = 'none'; % interpet alpha values as transparency values
        geoimg.FaceAlpha = 'texturemap'; % Indicate that the transparency can be different each pixel
        alpha(geoimg,double(~isnan(Topo_corrected)))
%         ax0.SortMethod = 'childorder';
        
          %Crameri color
%         cmap0=crameri('oleron');
%         colormap(ax0,cmap0);
%         clim([-4000 4000]);

          %ttcmap visualisation
%         [cmap,zlimits] = ttcmap(Topo_corrected,'cmap','mby');
%         colormap(ax0,cmap);
        
% Number of lines in the gradient
num_lines = 256;

% Interpolate additional values
% Topo_color = interp1(colors_topography.colors_topography(:, 1), colors_topography.colors_topography(:, 2:end), linspace(min(colors_topography.colors_topography(:, 1)), max(colors_topography.colors_topography(:, 1)), num_lines));
Topo_color = interp1(Topography_colors(:, 1),Topography_colors(:, 2:end), linspace(min(Topography_colors(:, 1)), max(Topography_colors(:, 1)), num_lines));
colormap(ax0,Topo_color./255); clim([-8000 4000]);

        c0=colorbar('Position',[0.1 0.18 0.02 0.66]);
        c0.Label.String= "Elevations [m]";
        c0.Label.FontSize = 14;
  
        hold on;
        ax2 = gca;
        setm(ax2, 'MapProjection', 'robinson')
         geoshow(Yeq, Xeq,continents_reversed_position_for_plot_dark,'facealpha', 0);
         contourm(Yeq,Xeq, continents, 0.3, 'k','LineWidth',2);

        contourm(Yeq,Xeq, continents_init, 0.3, 'w','LineWidth',2);
%         quivermc(Yeq_lr, Xeq_lr, V_theta_vec_surface, V_phi_vec_surface,'color','k','linewidth',1,'units','Velocity (m.yr^-^1)','reference',scaling_velocity_factor_surface);
%             plume_str = ['Plumes = ' num2str(numPlumes)];
%             subduction_str = ['Subductions = ' num2str(numSubductions)];
%             title([{time_str} plume_str subduction_str]);
            title({time_str});
        ax2.SortMethod='childorder';

%% Plot corrected topography
%  %% Plot Strain rate with topography
%  Topo_corrected=Topo-2200;
%  Topo_corrected(oceans) = ((Topo < 0) .* (1000 / 3300) .* abs(Topo));
% 
% %        figure();
%         figure(2);
%         set(gcf, 'Units', 'Normalized', 'OuterPosition', [0, 0.04, 1, 0.96],'color','w');clf;
%         set(gcf, 'color', 'w');
%         view(0, 90);
%         axesm mercator
%         framem;
%         gridm;
%         ax0 = gca;
%         setm(ax0, 'MapProjection', 'robinson')
%         geoshow(Yeq, Xeq,Topo_corrected,'facealpha', 0.3);
%         colormap(ax0);
%         ax0.SortMethod = 'childorder';
%         cmap0=crameri('oleron');
%         colormap(ax0,cmap0);
%         clim([-6000 6000]);
% 
% %         clim([-26000 -22000]);
%         % %caxis([-6000 6000])
%         % %         caxis([mean_surface_topography-abs(max_surface_topography-min_surface_topography)/2 mean_surface_topography+abs(max_surface_topography-min_surface_topography)/2])
% 
%         c0=colorbar('Position',[0.1 0.18 0.02 0.66]);
%         c0.Label.String= "Elevations [m]";
%         c0.Label.FontSize = 14;
%         freezeColors;
%         % c0=colorbar(ax0);
%         %%
%         hold on;
%         Plot continents with transparency
%         ax2 = gca;
%         setm(ax2, 'MapProjection', 'robinson')
%         geoshow(Yeq, Xeq,continents_reversed_position_for_plot,'facealpha', 0.3);
%         colormap(ax2);
%         contourm(Yeq,Xeq, continents_init, 0.3, 'w','LineWidth',2);
%         quivermc(Yeq_lr, Xeq_lr, V_theta_vec_surface, V_phi_vec_surface,'color','k','linewidth',1,'units','Velocity (m.yr^-^1)','reference',scaling_velocity_factor_surface);
%         ax2.SortMethod='childorder';

%% Track plumes and subductions over depths

        plumes_depths_tracking=[1,1,5];
        %select the depth files 1 is the shallowest and get deeper as the value increase
%         in the case where subduction_depths_tracking is 0, then the topography will be used to track the subductions
        subduction_depths_tracking=[0,1,5];
        plumes_non_adiabtic_tracking_temperature = [30,30,50];%5
        subduction_non_adiabtic_tracking_temperature = [0,-200,-400];%5
        trenches_elevation_threshold = -3000;

        % Create the histogram
%         figure();
%         bin_edges_topo = -8000:10:8000;
%         bin_counts = histcounts(Topo_corrected, bin_edges_topo);
%         histogram('BinEdges', bin_edges_topo, 'BinCounts', bin_counts, 'FaceColor', 'blue', 'EdgeColor', 'black');
%         xlabel('Topography');
%         ylabel('Frequency');
%         title('Histogram of Topography');

        
%         figure();pcolor(trenches);shading flat;

        for tt= 1:numel(plumes_depths_tracking)
         
            figure(1+tt);clf;
            set(gcf, 'Units', 'Normalized', 'OuterPosition', [0, 0.04, 1, 0.96],'color','w');clf;
            set(gcf, 'color', 'w');
            view(0, 90);
            ax1 = gca;
            axesm mercator
            framem;
            gridm;
            %         geoshow(Yeq, Xeq,non_adiabT_depths{displayed_depth}, 'DisplayType', 'surface');
            %         ax1 = gca;
            %         setm(ax1, 'MapProjection', 'robinson')
            %         %             zlabel('Surface');
            %         cmap1=crameri('vik');
            %         colormap(ax1,cmap1)
            %         c1 = colorbar();%ax1,'Position',[0.91 0.18 0.02 0.66]
            %         clim([-200 200]);
            %         c1.Label.String = ['Non adiabatic temperature [K] at ' num2str(unique_depths(displayed_depth)./1e3) ' km'];
            %         c1.Label.FontSize = 14;
            % set(ax1,'ColorScale');
            % % ticks = linspace((-5e-7), (5e-7), 6);
            % % tickLabels = arrayfun(@(x) sprintf('{%0.1f}', x), ticks, 'UniformOutput', false);
            % % c1.Ticks = ticks;
            % % c1.TickLabels = tickLabels;0

%             freezeColors(ax1);
            plumes_position = non_adiabT_depths{plumes_depths_tracking(tt)}>plumes_non_adiabtic_tracking_temperature(tt);
            plumes_positions2=1-plumes_position;
            plumes_position3 = ind2rgb(plumes_positions2 + 1, CMap_plume);
            geoshow(Yeq, Xeq, plumes_position3);
            setm(ax1, 'MapProjection', 'robinson')
            freezeColors(ax1);
            hold on;
             if (subduction_depths_tracking(tt)==0)
                  disp('here1');
                  current_depth = 0;
                  subductions_position = Topo_corrected<=trenches_elevation_threshold;
                  trenches_position2=1-subductions_position;
                  subductions_position3 = ind2rgb(trenches_position2 + 1, CMap_subduction)
                  cmatrix_subduction = contourm(Yeq,Xeq, Topo_corrected, [-90000 trenches_elevation_threshold], 'b','LineWidth',2);
             elseif (subduction_depths_tracking(tt)~=0)
                  disp('here2');
                  current_depth=unique_depths(subduction_depths_tracking(tt));
                  subductions_position = non_adiabT_depths{subduction_depths_tracking(tt)}<subduction_non_adiabtic_tracking_temperature(tt);
                  subductions_position2=1-subductions_position;
                  subductions_position3 = ind2rgb(subductions_position2 + 1, CMap_subduction);
                  cmatrix_subduction = contourm(Yeq,Xeq, non_adiabT_depths{subduction_depths_tracking(tt)}, [-90000 subduction_non_adiabtic_tracking_temperature(tt)], 'b','LineWidth',2);
             end            
            geoshow(Yeq, Xeq,subductions_position3,'FaceAlpha',0.5);
            boundaries_position3 = ind2rgb(Plate_boundaries2+1, CMap_boundaries);
            geoshow(Yeq, Xeq, boundaries_position3,'FaceAlpha',0.4);
%             setm(ax1, 'MapProjection', 'robinson')
            cmatrix_plume = contourm(Yeq,Xeq, non_adiabT_depths{plumes_depths_tracking(tt)}, [plumes_non_adiabtic_tracking_temperature(tt) plumes_non_adiabtic_tracking_temperature(tt)], 'r','LineWidth',2);
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
            else
                [x_subduction, y_subduction] = C2xyz(cmatrix_subduction);
                polyin_subductions = polyshape(x_subduction, y_subduction);
                numSubductions = polyin_subductions.NumRegions;
                It2 = bwmorph(subductions_position,'thin','inf');
%                 If some cleaning is needed
%                 It2 = bwmorph(It,'clean');
                Subduction_connectivity= bwconncomp(It2);
                numSubduction_connected = Subduction_connectivity.NumObjects;
               % This additional figure can display the endbranch of every subductions
%             figure();
%             imshow(It2);
%             [i_dist,j_dist] = find(bwmorph(It2,'endpoints'));
%             for n = 1:numel(i_dist)
%                 text(j_dist(n),i_dist(n),[num2str(D(i_dist(n),j_dist(n)))],'color','g');
%             end
            end
%             Give the length of each skeloteneize subduction in degree
%             test_length=regionprops(It2, 'area');

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


            %             To create plogygons
            %         poly_plume = struct('Geometry', 'Polygon', 'Longitude', x_plume, 'Latitude', y_plume);
            %                 poly_subduction = struct('Geometry', 'Polygon', 'Longitude', x_subduction, 'Latitude', y_subduction);
            %         quivermc(Yeq_lr, Xeq_lr, V_theta_vec_surface, V_phi_vec_surface,'color','k','linewidth',1,'units','Velocity (m.yr^-^1)','reference',scaling_velocity_factor_surface);
            %                 geoshape_plume=geoshape(poly_plume);
            %                  geoshape_subduction=geoshape(poly_subduction);
            %         legend_entry = sprintf('Plumes NonAdiabT > %d', plumes_non_adiabtic_tracking_temperature);
            %         legend_entry = sprintf('Subductions NonAdiabT < %d', subduction_non_adiabtic_tracking_temperature);
%             quivermc(Yeq_lr, Xeq_lr,V_theta_vec_depths{plumes_depths_tracking(tt)}, V_phi_vec_depths{plumes_depths_tracking(tt)},'color','k','linewidth',1,'units','Velocity (m.yr^-^1)','reference',scaling_velocity_factor_depths);
            plume_str = ['Plumes = ' num2str(numPlumes)];
            title([{time_str} plume_str]);
            subduction_str = ['Subductions = ' num2str(numSubductions)];
            title([{time_str} plume_str subduction_str]);
            contourm(Yeq,Xeq, continents_init, 0.3, 'k','LineWidth',2);
            ax2 = gca;
            setm(ax2, 'MapProjection', 'robinson')
            geoshow(Yeq, Xeq,continents_reversed_position_for_plot,'facealpha', 0.3);
            colormap(ax2);
            ax2.SortMethod = 'childorder';

            %%

            time_for_output_structures = surface.Time(1)./1e6;

            % Check if the table exists, if not, create a new table
            if ~exist('Geodynamics_structures_distribution', 'var')
                Geodynamics_structures_distribution = table(time_for_output_structures, numPlumes, numSubductions,numSubduction_connected,current_depth./1e3,Maxferet_mean,Maxferet_max,Maxferet_min,nferet, ...
                    'VariableNames', {'Time', 'NumPlumes', 'NumSubductions','NumSubduction_connected','Depths','Maxferet_mean','Maxferet_max','Maxferet_in','nferet'});
            else
                % Add new values to the existing table
                newRow = table(time_for_output_structures, numPlumes, numSubductions,numSubduction_connected,current_depth./1e3,Maxferet_mean,Maxferet_max,Maxferet_min,nferet, ...
                    'VariableNames', {'Time', 'NumPlumes', 'NumSubductions','NumSubduction_connected','Depths','Maxferet_mean','Maxferet_max','Maxferet_in','nferet'});

                Geodynamics_structures_distribution = vertcat(Geodynamics_structures_distribution, newRow);
            end

            %         % Display the updated table
            disp(Geodynamics_structures_distribution);

            % Save the updated table to a file (adjust the filename and format as needed)
            writetable(Geodynamics_structures_distribution, ['Distribution_' model_titles{m} '.csv']);
        end


        %% Plot Strain rate
        figure(2+tt);
%         clf;
        set(gcf, 'Units', 'Normalized', 'OuterPosition', [0, 0.04, 1, 0.96]);clf;
        set(gcf, 'color', 'w');
        view(0, 90);
        axesm mercator

        geoshow(Yeq, Xeq, log10(strain_rate), 'DisplayType', 'surface');
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
        ax2 = gca;
        setm(ax2, 'MapProjection', 'robinson')
        geoshow(Yeq, Xeq,continents_reversed_position_for_plot,'facealpha', 0.3);
        colormap(ax2);freezeColors;
        contourm(Yeq,Xeq, continents_init, 0.3, 'w','LineWidth',2);
        quivermc(Yeq_lr, Xeq_lr, V_theta_vec_surface, V_phi_vec_surface,'color','k','linewidth',1,'units','Velocity (m.yr^-^1)','reference',scaling_velocity_factor_surface);
        ax2.SortMethod='childorder';


        %%

        %         quiverm(Yeq,Xeq,V_theta,V_phi,'LineSpec',s=0,'filled');

        %         for cc = 1 : number_of_init_continent_boundaries
        %           thisBoundary = continent_init_boundaries{cc};
        %           contourf(thisBoundary(:,1)-91,thisBoundary(:,2)-181);
        %         end


        %% Plot Strain rate with topography
        % Get the active plate boundaries using a threshold on strain rate
        strain_rate_threshold = strain_rate;
        strain_rate_threshold(strain_rate_threshold<1e-15)=NaN;
        figure(3+tt);clf;
        %         set(gcf, 'Units', 'Normalized', 'OuterPosition', [0, 0.04, 1, 0.96],'color','w');clf;
        set(gcf, 'color', 'w');
        view(0, 90);
        axesm mercator
        framem;
        gridm;

        geoshow(Yeq,Xeq,Topo,'DisplayType', 'surface');
        ax0 = gca;
        setm(ax0, 'MapProjection', 'robinson')
        cmap0=crameri('oleron');
        colormap(ax0,cmap0);
        clim([-26000 -22000]);
        c0=colorbar('Position',[0.1 0.18 0.02 0.66]);
        c0.Label.String= "Elevations [m]";
        c0.Label.FontSize = 14;
        freezeColors;
        % c0=colorbar(ax0);

        hold on;

        geoshow(Yeq, Xeq, log10(strain_rate_threshold), 'DisplayType', 'surface');
        ax1 = gca;
        setm(ax1, 'MapProjection', 'robinson')
        zlabel('Surface');
        title({'Surface', time_str});
        cmap1=crameri('-roma');
        colormap(ax1,cmap1)
        c1 = colorbar(ax1,'Position',[0.91 0.18 0.02 0.66]);
        clim([-16 -14]);
        set(ax1,'ColorScale');
        ticks = linspace(log10(1e-16), log10(1e-14), 6);
        tickLabels = arrayfun(@(x) sprintf('{%0.1f}', x), ticks, 'UniformOutput', false);
        c1.Ticks = ticks;
        c1.TickLabels = tickLabels;
        c1.Label.String = 'Strain rate [1e^ /s]';
        c1.Label.FontSize = 14;
        freezeColors;


        % freezeColors(c1)
        % Plot continents with transparency
        ax2 = gca;
        setm(ax2, 'MapProjection', 'robinson')
        geoshow(Yeq, Xeq,continents_reversed_position_for_plot,'facealpha', 0.3);
        colormap(ax2);
        contourm(Yeq,Xeq, continents_init, 0.3, 'w','LineWidth',2);
        quivermc(Yeq_lr, Xeq_lr, V_theta_vec_surface, V_phi_vec_surface,'color','k','linewidth',1,'units','Velocity (m.yr^-^1)','reference',scaling_velocity_factor_surface);
        ax2.SortMethod='childorder';
        
        %% Plot oceanic age 
%         idx_oceanic_domain = find(oceanic == 0);
        oceanic_age_domain = oceanic_age;
        oceanic_age_domain(continents_threshold) = NaN;
        % Get the active plate boundaries using a threshold on strain rate
        figure(4+tt);
%         clf;
        set(gcf, 'Units', 'Normalized', 'OuterPosition', [0, 0.04, 1, 0.96],'color','w');clf;
        set(gcf, 'color', 'w');
        view(0, 90);
        axesm mercator
        framem;
        gridm;

        geoshow(Yeq, Xeq, oceanic_age_domain, 'DisplayType', 'surface');
        ax1 = gca;
        setm(ax1, 'MapProjection', 'robinson')
        zlabel('Surface');
        title({'Surface', time_str});
        cmap1=crameri('-roma');
        colormap(ax1,cmap1)
        c1 = colorbar(ax1,'Position',[0.91 0.18 0.02 0.66]);
        clim([0 200]);
        set(ax1,'ColorScale');
        ticks = linspace(0, 200, 6);
        tickLabels = arrayfun(@(x) sprintf('{%0.1f}', x), ticks, 'UniformOutput', false);
        c1.Ticks = ticks;
        c1.TickLabels = tickLabels;
        c1.Label.String = 'Oceanice age [My]';
        c1.Label.FontSize = 14;
        freezeColors;


        % freezeColors(c1)
        % Plot continents with transparency
        ax2 = gca;
        setm(ax2, 'MapProjection', 'robinson')
        geoshow(Yeq, Xeq, continents_reversed_position_for_plot,'facealpha', 0.3);
        colormap(ax2);
        contourm(Yeq,Xeq, continents_init, 0.3, 'w','LineWidth',2);
        quivermc(Yeq_lr, Xeq_lr, V_theta_vec_surface, V_phi_vec_surface,'color','k','linewidth',1,'units','Velocity (m.yr^-^1)','reference',scaling_velocity_factor_surface);
        ax2.SortMethod='childorder';

        
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


        %% Plot divergence
        smooth_divergence = imgaussfilt(divergence,0.2);
        figure(5+tt);
%         clf;
        set(gcf, 'Units', 'Normalized', 'OuterPosition', [0, 0.04, 1, 0.96],'color','w');
        clf;set(gcf, 'color', 'w');
        view(0, 90);
        axesm mercator
        geoshow(Yeq, Xeq,smooth_divergence, 'DisplayType', 'surface');
        ax1 = gca;
        setm(ax1, 'MapProjection', 'robinson')
        zlabel('Surface');
        title({'Surface', time_str});
        cmap1=crameri('vik');
        colormap(ax1,cmap1)
        c1 = colorbar();%ax1,'Position',[0.91 0.18 0.02 0.66]
        clim([-5e-8 5e-8]);
        set(ax1,'ColorScale');
        % ticks = linspace((-5e-7), (5e-7), 6);
        % tickLabels = arrayfun(@(x) sprintf('{%0.1f}', x), ticks, 'UniformOutput', false);
        % c1.Ticks = ticks;
        % c1.TickLabels = tickLabels;
        c1.Label.String = 'Divergence [1e^ /s]';
        c1.Label.FontSize = 14;
        framem;
        gridm;
        hold on;
        % Plot continents with transparency
        geoshow(Yeq, Xeq, continents_reversed_position_for_plot,'facealpha', 0.3);
        colormap();
        contourm(Yeq,Xeq, continents_init, 0.3, 'w','LineWidth',2);
        quivermc(Yeq_lr, Xeq_lr, V_theta_vec_surface, V_phi_vec_surface,'color','k','linewidth',1,'units','Velocity (m.yr^-^1)','reference',scaling_velocity_factor_surface);
        SortMethod='childorder';
        
        %% Plot Adibiatic temperature at different depths

        for uuu = 1:numel(unique_depths)
            figure(5+tt+uuu);
%             clf;
            set(gcf, 'Units', 'Normalized', 'OuterPosition', [0, 0.04, 1, 0.96],'color','w');clf;
            set(gcf, 'color', 'w');
            view(0, 90);
            axesm mercator
            geoshow(Yeq, Xeq,non_adiabT_depths{uuu}, 'DisplayType', 'surface');
            ax1 = gca;
            setm(ax1, 'MapProjection', 'robinson')
            %             zlabel('Surface');
            title({time_str});
            cmap1=crameri('vik');
            colormap(ax1,cmap1)
            c1 = colorbar();%ax1,'Position',[0.91 0.18 0.02 0.66]
            clim([-200 200]);
            % set(ax1,'ColorScale');
            % % ticks = linspace((-5e-7), (5e-7), 6);
            % % tickLabels = arrayfun(@(x) sprintf('{%0.1f}', x), ticks, 'UniformOutput', false);
            % % c1.Ticks = ticks;
            % % c1.TickLabels = tickLabels;0
            c1.Label.String = ['Non adiabatic temperature [K] at ' num2str(unique_depths(uuu)./1e3) ' km'];
            c1.Label.FontSize = 14;
            framem;
            gridm;
            hold on;
            ax2 = gca;
            setm(ax2, 'MapProjection', 'robinson')
            geoshow(Yeq, Xeq,continents_reversed_position_for_plot,'facealpha', 0.3);
            colormap(ax2);freezeColors;
            contourm(Yeq,Xeq, continents_init, 0.3, 'w','LineWidth',2);
            quivermc(Yeq_lr, Xeq_lr,V_theta_vec_depths{uuu}, V_phi_vec_depths{uuu},'color','k','linewidth',1,'units','Velocity (m.yr^-^1)','reference',scaling_velocity_factor_depths);
            ax2.SortMethod = 'childorder';
        end
        %%
        for vvv = 1:numel(unique_depths)
            %%
            figure(5+tt+uuu+vvv);
%             clf;
            set(gcf, 'Units', 'Normalized', 'OuterPosition', [0, 0.04, 1, 0.96],'color','w');clf;
            set(gcf, 'color', 'w');
            view(0, 90);
            axesm mercator
            geoshow(Yeq, Xeq,V_magnitude_depths{vvv}, 'DisplayType', 'surface');
            ax1 = gca;
            setm(ax1, 'MapProjection', 'robinson')
            %             zlabel('Surface');
            title({time_str});
            cmap1=crameri('vik');
            colormap(ax1,cmap1)
            c1 = colorbar();%ax1,'Position',[0.91 0.18 0.02 0.66]
            clim([0 0.05]);
            % set(ax1,'ColorScale');
            % % ticks = linspace((-5e-7), (5e-7), 6);
            % % tickLabels = arrayfun(@(x) sprintf('{%0.1f}', x), ticks, 'UniformOutput', false);
            % % c1.Ticks = ticks;
            % % c1.TickLabels = tickLabels;0
            c1.Label.String = ['Velocity magnitude [m.yr^-^1] at ' num2str(unique_depths(vvv)./1e3) ' km'];
            c1.Label.FontSize = 14;
            framem;
            gridm;
            hold on;
            ax2 = gca;
            setm(ax2, 'MapProjection', 'robinson')
            geoshow(Yeq, Xeq,continents_reversed_position_for_plot,'facealpha', 0.3);
            colormap(ax2);freezeColors;
            contourm(Yeq,Xeq, continents_init, 0.3, 'w','LineWidth',2);
            quivermc(Yeq_lr, Xeq_lr,V_theta_vec_depths{vvv}, V_phi_vec_depths{vvv},'color','k','linewidth',1,'units','Velocity (m.yr^-^1)','reference',scaling_velocity_factor_depths);
            ax2.SortMethod='childorder';

            %%
        end

        %         Calculate regional Mobility = lithosphere vel / mantle
              index_depth_vel_rms=3 ; %         depth vel at 1000km
        Regional_Mobility = V_magnitude_surface./(V_magnitude_depths{3});

        figure(6+tt+uuu+vvv);clf;
        set(gcf, 'Units', 'Normalized', 'OuterPosition', [0, 0.04, 1, 0.96],'color','w');clf;
        set(gcf, 'color', 'w');
        view(0, 90);
        axesm mercator
        geoshow(Yeq, Xeq,Regional_Mobility, 'DisplayType', 'surface');
        ax1 = gca;
        setm(ax1, 'MapProjection', 'robinson')
        %             zlabel('Surface');
        title({time_str});
        cmap1=crameri('vik');
        colormap(ax1,cmap1)
        c1 = colorbar();%ax1,'Position',[0.91 0.18 0.02 0.66]
        clim([0 5]);
        % set(ax1,'ColorScale');
        % % ticks = linspace((-5e-7), (5e-7), 6);
        % % tickLabels = arrayfun(@(x) sprintf('{%0.1f}', x), ticks, 'UniformOutput', false);
        % % c1.Ticks = ticks;
        % % c1.TickLabels = tickLabels;0
        c1.Label.String = ['Regional mobility V_surface / V_' num2str(unique_depths(index_depth_vel_rms)./1e3) ' km'];
        c1.Label.FontSize = 14;
        framem;
        gridm;
        hold on;
        ax2 = gca;
        setm(ax2, 'MapProjection', 'robinson')
        geoshow(Yeq, Xeq,continents_reversed_position_for_plot,'facealpha', 0.3);
        colormap(ax2);freezeColors;
        contourm(Yeq,Xeq, continents_init, 0.3, 'w','LineWidth',2);
        quivermc(Yeq_lr, Xeq_lr, V_theta_vec_surface, V_phi_vec_surface,'color','k','linewidth',1,'units','Velocity (m.yr^-^1)','reference',scaling_velocity_factor_surface);
        ax2.SortMethod='childorder';



        %% Video in series (in construction)
        if strcmp(make_video, 'true')
            % Capture the current frame
            frame = getframe(gcf);

            % Write the frame to the video file
            writeVideo(videoWriter, frame);
        else 
            tilefigs;
        end

    end

    if strcmp(make_video, 'true')
        % Close the video file
        close(videoWriter);
    end
end