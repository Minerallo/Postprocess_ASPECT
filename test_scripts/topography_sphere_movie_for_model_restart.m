clc ; close all; clear all;
% '/Users/ponsm/Desktop/modelblogin/model/globalscale/sphere3d/P01a_Pangea_1GPa_Mantle_C40MPa_LR/',
% '/Users/ponsm/Desktop/modelblogin/model/globalscale/sphere3d/P01b_Pangea_1GPa_Mantle_C40MPa_HR/',
path_models = {'/Users/ponsm/Desktop/modelblogin/model/globalscale/sphere3d/P01c_Pangea_1GPa_Mantle_C60MPa_LR/'}

% Extracting model titles following 'sphere3d/'
model_titles = cell(size(path_models));
for t = 1:numel(path_models)
    splitted = strsplit(path_models{t}, '/');
    idx = find(contains(splitted, 'sphere3d'));
    if numel(idx) > 0
        model_titles{t} = splitted{idx+1};
    end
end


for m = 1:numel(path_models)
    path_model = cell2mat(path_models(m));
    make_video='false';
    resample_dynamic_topography = 10;

    statistic_parameters = {'Time'};
    % Get the data and corresponding statistic numbers and indices
    [data_stats, stats_number, header, stat_indices] = get_statistics(path_model, statistic_parameters);

    % Statistics
    time_vec = data_stats(:, 2)./1e6;
    time=num2str(time_vec, '%05.2f Myr');

    % find how many zeros for the initial refinement steps
    numbers_of_IRS=nnz(~time_vec);

    % Determine model length from first topography file
    filename = fullfile(path_model, 'topography.00000');
    dynamic_topography = import_dynamic_topography3D(filename);
    % modellength = max(dynamic_topography.x) / 1e3;

    file = dir(fullfile(path_model, 'topography.*'));
    filenames = {file.name};
    num_files = numel(filenames);

    % Specify the desired frame size for the video
    frameSize = [2000, 1000];

    % Create the figure with the specified size
    figure('WindowState', 'maximized','Position', [100 100 frameSize]);
    scatterSize =70;

    if strcmp(make_video, 'true')
        % Specify the video filename
        videoFilename = [model_titles{m} '.mp4'];

        % Create a VideoWriter object
        videoWriter = VideoWriter(videoFilename, 'MPEG-4');
        videoWriter.FrameRate = 4; % Specify the frame rate (adjust as needed)
        open(videoWriter);
    end

    for i = 2381:resample_dynamic_topography:num_files
        % for i =num_files
        % for i=
        fprintf('Processing dynamic topography file %d of %d\n', i, num_files)
        file_index = i;

        filename_extract = filenames{file_index};
        parts = split(filename_extract, '.');
        number_step = str2double(parts{2});
        time_str = sprintf('Time: %05.2f Myr', time_vec(number_step+numbers_of_IRS));

        % Import the data
        dynamic_topography = import_dynamic_topography3D(fullfile(path_model, filenames{file_index}));

        % Sort the data in the desired order
        x = dynamic_topography.x;
        y = dynamic_topography.y;
        z = dynamic_topography.z;
        surface_topography = dynamic_topography.surface_topography;

        % Calculate the radial distance from the center of the Earth
        r = sqrt(x.^2 + y.^2 + z.^2);

        % Calculate the azimuthal angle
        theta = atan2(y, x);

        % Sort the data based on the azimuthal angle
        [theta_sorted, theta_index] = sort(theta);
        x_sorted = x(theta_index);
        y_sorted = y(theta_index);
        z_sorted = z(theta_index);
        r_sorted = r(theta_index);

        surface_topography_sorted = surface_topography(theta_index);
        final_surface_topography_sorted = surface_topography_sorted - ((surface_topography_sorted < 2400) .* (1000 / 3300) .* abs(surface_topography_sorted - 2400)) - 2400;

        % Define the downsampling factor
        downsample_factor = 1;

        x_downsampled = x_sorted(1:downsample_factor:end);
        y_downsampled = y_sorted(1:downsample_factor:end);
        z_downsampled = z_sorted(1:downsample_factor:end);
        surface_topography_downsampled = final_surface_topography_sorted(1:downsample_factor:end);

        min_surface_topography = min(surface_topography_downsampled);
        max_surface_topography = max(surface_topography_downsampled);
        mean_surface_topography = mean(surface_topography_downsampled);

        % Calculate the extruded Z-coordinates based on the topography values
        extrusion_factor = 1; % Modify this value as desired
        z_extruded = z_downsampled * extrusion_factor;

        %     % Create a scatter plot with varying sizes and colors
        %     h1=subplot(1, 2, 1);
        %     scatter3(x_downsampled, y_downsampled, z_extruded, scatterSize, surface_topography_downsampled, 'filled');
        %     colorbar;
        %     xlabel('X');
        %     ylabel('Y');
        %     zlabel('Z');
        %     time_str = sprintf('Time: %05.2f Myr', time_vec(file_index));
        %     title(h1, {'Surface Topography', time_str});
        %     crameri('oleron');caxis([min_surface_topography min_surface_topography+20000]);%caxis([-6000 6000]);
        %     c=colorbar;set(gcf,'color','w');axis equal;c.Label.String= "Elevations [m]";
        %     axis off;  % Remove the axis lines
        %     grid off;  % Remove the grid lines

        % Rotate the view of subplot(1, 2, 1) at each iteration
        % rotation_angle_increment
        %     angle = (file_index - 1) * 5;
        %     view(h1, angle, 5);
        %  set(h1, 'CameraUpVector', [0 0 1]); % Fix the axes
        %

        % Convert Cartesian coordinates to longitude and latitude
        longitude_downsampled = atan2(y_downsampled, x_downsampled) * 180 / pi;
        latitude_downsampled = asin(z_downsampled ./ sqrt(x_downsampled.^2 + y_downsampled.^2 + z_downsampled.^2)) * 180 / pi;



        % Plot the data on a surface
        %     h2=subplot(1, 2, 2);
        % Create a grid of longitude and latitude values
        [Xeq, Yeq] = meshgrid(-180:1:180, -90:1:90);

        % Reshape the surface topography data to match the grid size
        Topo = griddata(longitude_downsampled, latitude_downsampled, surface_topography_downsampled, Xeq, Yeq);
        Topo = Topo(1:180,1:360);
        Xeq = Xeq(1:180,1:360);
        Yeq = Yeq(1:180,1:360);

        %     mapshow(Xeq,Yeq,Topo,"DisplayType","surface")
        %     surf(Xeq, Yeq, Topo, 'EdgeColor', 'none');
        %     colorbar;
        %     xlabel('Longitude');
        %     ylabel('Latitude');
        %     zlabel('Surface Topography');
        %     title(h2, {'Surface Topography', time_str});c=colorbar;crameri('oleron');caxis([min_surface_topography min_surface_topography+15000]);set(gcf,'color','w'); view(0,90);c.Label.String= "Elevations [m]";


        clf;
        set(gcf, 'Units', 'Normalized', 'OuterPosition', [0, 0.04, 1, 0.96]);
        axesm mercator
        geoshow(Yeq,Xeq,Topo,'DisplayType', 'surface')
        ax = gca;
        setm(ax,"MapProjection","robinson")
        colorbar;
        %     xlabel('Longitude');
        %     ylabel('Latitude');
        zlabel('Surface Topography');
        title({'Surface Topography', time_str});c=colorbar;crameri('oleron');set(gcf,'color','w'); view(0,90);c.Label.String= "Elevations [m]";
        caxis([mean_surface_topography-abs(max_surface_topography-min_surface_topography)/2 mean_surface_topography+abs(max_surface_topography-min_surface_topography)/2])
        framem
        gridm

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