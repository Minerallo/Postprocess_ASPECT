clc ; close all; clear all;
addpath(genpath('./Individual_scripts/'))

% path_models = {'/Users/ponsm/Desktop/modelblogin/model/globalscale/Potsprocess/sphere3d/P01a_Pangea_1GPa_Mantle_C40MPa_LR/'}
path_models = {'/Users/ponsm/Desktop/Model_sphere/Models_HLRN/P01a_Pangea_1GPa_Mantle_C40MPa_LR/'}
%     '/Users/ponsm/Desktop/Model_sphere/Models_HLRN/R01f_Rodinia_2GPa_Mantle_C10MPa_f005_LR/'}
%     '/Users/ponsm/Desktop/Model_sphere/Models_HLRN/P01a_Pangea_1GPa_Mantle_C40MPa_LR/'}

% Extracting model titles
model_titles = cell(size(path_models));
for t = 1:numel(path_models)
% Split the path into parts
parts = strsplit(path_models{t}, '/');
% Extract the last part
model_titles{t} = parts{end-1};
end

%Load colors
% colors_topography = load('colors.mat');

for m = 1:numel(path_models)
    path_model = cell2mat(path_models(m));
    make_video='make';
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
%     for i =99
            for i = 1:num_files_surface

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

%         theta_lithosphere = atan2(y_lithosphere, x_lithosphere);
%         [theta_sorted_lithosphere, theta_index_lithosphere] = sort(theta_lithosphere);
%         x_sorted_lithosphere = x_lithosphere(theta_index_lithosphere);
%         y_sorted_lithosphere = y_lithosphere(theta_index_lithosphere);
%         z_sorted_lithosphere = z_lithosphere(theta_index_lithosphere);
%         divergence_lithosphere_sorted = lithosphere.Divergence(theta_index_lithosphere);
%         oceanic_age_lithosphere_sorted = lithosphere.oceanic_age(theta_index_lithosphere);
% 
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

% 
%         x_surface_downsampled = x_sorted_surface(1:downsample_factor:end);
%         y_surface_downsampled = y_sorted_surface(1:downsample_factor:end);
%         z_surface_downsampled = z_sorted_surface(1:downsample_factor:end);
%         surface_topography_surface_downsampled = final_surface_topography_sorted(1:downsample_factor:end);

%         x_lithosphere_downsampled = x_sorted_lithosphere(1:downsample_factor:end);
%         y_lithosphere_downsampled = y_sorted_lithosphere(1:downsample_factor:end);
%         z_lithosphere_downsampled = z_sorted_lithosphere(1:downsample_factor:end);


        min_surface_topography = min(final_surface_topography_sorted);
        max_surface_topography = max(final_surface_topography_sorted);
        mean_surface_topography(count_ite) = mean(final_surface_topography_sorted);

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


%         % Convert Cartesian coordinates to longitude and latitude
        longitude_surface_downsampled = atan2(y_sorted_surface, x_sorted_surface) * 180 / pi;
        latitude_surface_downsampled = asin(z_sorted_surface ./ sqrt(x_sorted_surface.^2 + y_sorted_surface.^2 + z_sorted_surface.^2)) * 180 / pi;

%         longitude_lithosphere_downsampled = atan2(y_lithosphere_downsampled, x_lithosphere_downsampled) * 180 / pi;
%         latitude_lithosphere_downsampled = asin(z_lithosphere_downsampled ./ sqrt(x_lithosphere_downsampled.^2 + y_lithosphere_downsampled.^2 + z_lithosphere_downsampled.^2)) * 180 / pi;
% 


        % Plot the data on a surface
        %     h2=subplot(1, 2, 2);
        % Create a grid of longitude and latitude values
%         [Xeq, Yeq] = meshgrid(-180:1:180, -90:1:90);
        [Xeq_lr, Yeq_lr] = meshgrid(-180:5:180, -90:5:90);
        % Reshape the surface topography data to match the grid size
        Topo = griddata(longitude_surface_downsampled, latitude_surface_downsampled, final_surface_topography_sorted, Xeq, Yeq);
        Topo = Topo+diff_mean_surface_topography(count_ite);
%         V_phi_vec_surface = griddata(longitude_surface_downsampled, latitude_surface_downsampled,surface_v_phi_sorted, Xeq_lr, Yeq_lr);
%         V_theta_vec_surface = griddata(longitude_surface_downsampled, latitude_surface_downsampled,surface_v_theta_sorted, Xeq_lr, Yeq_lr);
%         V_phi_surface = griddata(longitude_surface_downsampled, latitude_surface_downsampled,surface_v_phi_sorted, Xeq, Yeq);
%         V_theta_surface = griddata(longitude_surface_downsampled, latitude_surface_downsampled,surface_v_theta_sorted, Xeq, Yeq);
%         V_r_surface = griddata(longitude_surface_downsampled, latitude_surface_downsampled,surface_v_r_sorted, Xeq, Yeq);
%         V_magnitude_surface = sqrt(V_r_surface.^2+V_phi_surface.^2+V_theta_surface.^2);
% 

%         strain_rate = griddata(longitude_surface_downsampled, latitude_surface_downsampled, surface_strain_rate_sorted, Xeq, Yeq);
        continents = griddata(longitude_surface_downsampled, latitude_surface_downsampled, surface_continent_sorted, Xeq, Yeq);

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

        
%% Plot corrected topography
 %% Plot Strain rate with topography
%  Topo_corrected=Topo;
 smooth_Topo = imgaussfilt(Topo,0.3);

 Topo_corrected=smooth_Topo ; 
 idx_oceans=find(oceans == 1);
  idx_cont=find(oceans == 0);
Topo_corrected(idx_oceans) = Topo_corrected(idx_oceans) - (Topo_corrected(idx_oceans)<0).* (1000 / 3300).* abs(Topo_corrected(idx_oceans));

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
% clim([-8000 4000]);
        
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
       
        
        %%
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