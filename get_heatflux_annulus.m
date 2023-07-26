function [time_heatflux, heatfluxmap, x_axis_interp_heatflux] = get_heatflux_annulus(path_model, dt_heatflux, resample_heatflux)

% Resample the model time step using the resample_topography variable
dt = dt_heatflux;
dt = dt * resample_heatflux;


% Determine model length from first topography file
filename = fullfile(path_model, 'heat_flux.00000');
topography = import_heatflux(filename);
modellength = max(topography.x) / 1e3;

file = dir(fullfile(path_model, 'heat_flux.*'));
filenames = {file.name};
num_files = numel(filenames);

% Create a new vector for interpolation with high resolution (10000 points)
resolution_heatflux=10000;
x_axis_interp_heatflux = 0:360/resolution_heatflux:360;



% Resample topo if needed for faster computing

% elevation_map = zeros(num_files, 10000 / (resample_topography * 1e3));
for i = 1:resample_heatflux:num_files
    fprintf('Processing heatflux file %d of %d\n', i, num_files)
    file_index = i;
%     if i == 1
%         file_index = 2;
%     end
    
    % Import the data
    topography = import_heatflux(fullfile(path_model, filenames{file_index}));

    % Sort the data in the desired order
    x = topography.VarName1;
    y = topography.x;
    z = topography.y;
    r = sqrt(x.^2 + y.^2);
    %theta allows to know in which quadrant of the circle we are. 
    theta = atan2(y, x); 
    %now we can sort the data depending on theta index
    [theta_sorted, theta_index] = sort(theta);
    x_sorted = x(theta_index);
    y_sorted = y(theta_index);
    r_sorted = r(theta_index);
    z_sorted = z(theta_index);
    
    % Convert theta to degrees and adjust the range to be 0 to 360 degrees
    theta_degrees = theta_sorted * 180 / pi;
    theta_degrees(theta_degrees<0) = theta_degrees(theta_degrees<0) + 360;
    
    % Find unique theta values and corresponding indices
    [theta_unique, index_unique] = unique(theta_degrees);
    
    
    % Interpolate the sorted z data to match the length of theta_degrees
    z_interp = interp1(theta_unique, z_sorted(index_unique), x_axis_interp_heatflux);
    elevation_map(i, :) = z_interp ;

end
heatfluxmap = elevation_map(:, :);
% time_heatflux=0:(dt*(num_files/resample_heatflux))/(num_files/resample_heatflux):dt*(num_files/resample_heatflux);
time_heatflux = 0:dt:dt * (size(heatfluxmap, 1) - 1);
end

%%%%check if topography is right 
% figure();
% plot(x_sorted, y_sorted, 'LineWidth', 2);
% figure();
% plot(theta_degrees, r_sorted-6371e3, 'LineWidth', 2);
% hold on;
% plot(theta_degrees, z_sorted, 'LineWidth', 2);
% plot(vec_interp, z_interp, 'LineWidth', 2);
% xlabel('Angle (degrees)');
% ylabel('Topography (m)');
% legend('Topography recalculated from coordinates (m)', 'Aspect topography (m)');
% figure();
% 
% plot(vec_interp, z_interp, 'LineWidth', 2);
% xlabel('Angle (degrees)');
% ylabel('Topography (m)');
% legend('Aspect topography (m)');




% %% %%%%%%%%%%
% if any(Output_numbers==[1,3,8])
%     %Plot every 10 Ma
%     dt_elevation=10e6/dt;
%     
%     figure(3); clf;
%     % plot(topogreatcircle(:,3)-topogreatcircle(indexmin_topo,3),topogreatcircle(:,end));
%     plot(topogreatcircle(:,3)-topogreatcircle(indexmin_topo,3),topocorrected/1e3,'LineWidth',line_width);hold on ;
%     for i=1:fix(size(Elevationmap,1)/dt_elevation)
%         plot(Vec_lenght-Vec_lenght(indexmin_model(1)),(Elevationmap(i*dt_elevation,:)-bathymetry)/1e3,'LineWidth',line_width);
%         
%     end
%     xlabel('Distance from the trench [km]');
%     ylabel('Altitude [km]'); title('Evolution of surface topography'); set(gcf,'color','w');%set(gca, 'color', 'none');
%     legend ('Topography corrected from water','10 Ma','20 Ma','30 Ma','40 Ma','50 Ma');
% end

% figure(3);clf;
% plot(topogreatcircle(:,3)-topogreatcircle(indexmin_topo,3),topocorrected/1e3,'LineWidth',line_width);hold on ;
% plot(Vec_lenght-Vec_lenght(indexmin_model(1)),(Elevationmap(end,:)-bathymetry)/1e3,'LineWidth',line_width);
% legend('Real topography corrected from water','Model');set(gcf,'color','w');
% xlabel('Distance from the trench [km]');