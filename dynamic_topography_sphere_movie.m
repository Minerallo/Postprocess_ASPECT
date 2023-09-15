clc ; close all; clear all;
% path_model ='/Users/ponsm/Desktop/modelblogin/model/globalscale/sphere3d/3D_ascii/013_012_initglobal_restart/';
% path_model ='/Users/ponsm/Desktop/bbk00014/globalscale/sphere/019_013_init_continents_density_continents_only_maxtinestp_faster/'
% path_model ='/Users/ponsm/Desktop/bbk00014/globalscale/sphere/021c_013_cont_only_plate/'
% path_model ='/Users/ponsm/Desktop/bbk00014/globalscale/sphere/023_test_vel_and_restart/';
path_model ='/Users/ponsm/Desktop/modelblogin/model/globalscale/sphere3d/030_geometric_gplate_fast_1e21_kinematic';
make_video='false';
resample_dynamic_topography = 5;

statistic_parameters = {'Time'};
% Get the data and corresponding statistic numbers and indices
[data_stats, stats_number, header, stat_indices] = get_statistics(path_model, statistic_parameters);

% Statistics
time_vec = data_stats(:, 2)./1e6;
time=num2str(time_vec, '%05.2f Myr');

% Determine model length from first topography file
filename = fullfile(path_model, 'dynamic_topography_surface.00000');
dynamic_topography = import_dynamic_topography3D(filename);
modellength = max(dynamic_topography.x) / 1e3;

file = dir(fullfile(path_model, 'dynamic_topography_surface.*'));
filenames = {file.name};
num_files = numel(filenames);

% Specify the desired frame size for the video
frameSize = [2000, 1000];

% Create the figure with the specified size
figure('WindowState', 'maximized','Position', [100 100 frameSize]);
scatterSize =70;

if strcmp(make_video, 'true')
    % Specify the video filename
    videoFilename = 'topography_rotation_video.mp4';

    % Create a VideoWriter object
    videoWriter = VideoWriter(videoFilename, 'MPEG-4');
    videoWriter.FrameRate = 4; % Specify the frame rate (adjust as needed)
    open(videoWriter);
end

for i = num_files
%     1:resample_dynamic_topography:num_files
    fprintf('Processing dynamic topography file %d of %d\n', i, num_files)
    file_index = i;

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

    % Calculate the extruded Z-coordinates based on the topography values
    extrusion_factor = 1; % Modify this value as desired
    z_extruded = z_downsampled * extrusion_factor;

    % Create a scatter plot with varying sizes and colors
    h1=subplot(1, 2, 1);
    scatter3(x_downsampled, y_downsampled, z_extruded, scatterSize, surface_topography_downsampled, 'filled');
    colorbar;
    xlabel('X');
    ylabel('Y');
    zlabel('Z');
    time_str = sprintf('Time: %05.2f Myr', time_vec(file_index));
    title(h1, {'Surface Topography', time_str});
    crameri('oleron');caxis([-6000 6000]);c=colorbar;set(gcf,'color','w');axis equal;c.Label.String= "Elevations [m]";

    % Rotate the view of subplot(1, 2, 1) at each iteration
    % rotation_angle_increment
%     angle = (file_index - 1) * 5;
%     view(h1, angle, 5);
%  set(h1, 'CameraUpVector', [0 0 1]); % Fix the axes
    axis off;  % Remove the axis lines
    grid off;  % Remove the grid lines

    % Convert Cartesian coordinates to longitude and latitude
    longitude_downsampled = atan2(y_downsampled, x_downsampled) * 180 / pi;
    latitude_downsampled = asin(z_downsampled ./ sqrt(x_downsampled.^2 + y_downsampled.^2 + z_downsampled.^2)) * 180 / pi;



    % Plot the data on a surface
    h2=subplot(1, 2, 2);
    % Create a grid of longitude and latitude values
    [Xeq, Yeq] = meshgrid(-180:1:180, -90:1:90);

    % Reshape the surface topography data to match the grid size
    Topo = griddata(longitude_downsampled, latitude_downsampled, surface_topography_downsampled, Xeq, Yeq);

    surf(Xeq, Yeq, Topo, 'EdgeColor', 'none');
    colorbar;
    xlabel('Longitude');
    ylabel('Latitude');
    zlabel('Surface Topography');
    title(h2, {'Surface Topography', time_str});c=colorbar;crameri('oleron');caxis([-6000 6000]);set(gcf,'color','w'); view(0,90);c.Label.String= "Elevations [m]";
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