function [time_elevation, elevation, x_axis_interp, dip_topography] = get_topography_annulus(path_model, dt_topography, resample_topography, calculate_topography_dip,topography_smoothing_interval_for_dip_calculation)

% Resample the model time step using the resample_topography variable
dt = dt_topography;
dt = dt * resample_topography;


% Determine model length from first topography file
filename = fullfile(path_model, 'topography.00000');
topography = import_topography(filename);
modellength = max(topography.x) / 1e3;

file = dir(fullfile(path_model, 'topography.*'));
filenames = {file.name};
num_files = numel(filenames);

% Create a new vector for interpolation with high resolution (10000 points)
resolution_topo=10000;
x_axis_interp = 0:360/resolution_topo:360;



% Resample topo if needed for faster computing
count_ite=0;
% elevation_map = zeros(num_files, 10000 / (resample_topography * 1e3));
for i = 1:resample_topography:num_files
    count_ite=count_ite+1;

    % Check if the next iteration will go over num_files
    if (i > num_files)
        break;  % Exit the loop if the next iteration exceeds num_files
    end

    fprintf('Processing topography file %d of %d\n', i, num_files)
    file_index = i;
%     if i == 1
%         file_index = 2;
%     end
    
    % Import the data
    topography = import_topography(fullfile(path_model, filenames{file_index}));

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
    z_interp = interp1(theta_unique, z_sorted(index_unique), x_axis_interp);
    elevation_map(count_ite, :) = z_interp ;

    %     if dip calculation is asked
    if strcmp(calculate_topography_dip, 'true')

    % Check if the next iteration will go over num_files
    if i > num_files
        break;  % Exit the loop if the next iteration exceeds num_files
    end

%    find the indices of the places where there are no elevations2
%    (depracted from topolayer might not be relevant here)
            index_nan=find(isnan(z_interp));
%    The topography might be noisy because of the adaptative mesh refinement 
%    we can apply a slight smoothing to it
            smooth_topo=smooth(z_interp,topography_smoothing_interval_for_dip_calculation)';
%    Replace nan value by zero after smoothing (depracted from topolayer might not be relevant here)
            smooth_topo(index_nan)=0;
            dipsmooth=atand((diff(smooth_topo)./1e3));
            dip_map(count_ite,:)=dipsmooth;
%    For comparison with no smoothing uncomment the following lines
            z_interp(index_nan)=0;
            diptopo=atand((diff(z_interp)./1e3));
%    The figure will plot the last profil for comparison with and without smoothing
            figure();
            plot(1:numel(diptopo), diptopo, 'r-', 'LineWidth', 1.4);
            hold on;
            plot(1:numel(dipsmooth), dipsmooth, 'b-', 'LineWidth', 1.4);
            legend('Dip Not Smoothed','Dip Smoothed');
            xlabel('Model length');
            ylabel('Dip °');title('Dip Topography Comparison');legend('Location', 'Best');  
            
    end  

end
elevation = elevation_map(:, :);
time_elevation = 0:dt:dt * (size(elevation, 1) - 1);
% time_elevation=0:(dt*(num_files/resample_topography))/(num_files/resample_topography):dt*(num_files/resample_topography); %not working 

if strcmp(calculate_topography_dip, 'true')
    dip_topography = dip_map(:, :);
else
    dip_topography = 0;
end

end

%%%%% A new module in construction so the iteration ver the file does not
%%%%% go beyound the number of file when resampling
% function [time_elevation, elevation, x_axis_interp] = get_topography_annulus(path_model, dt_topography, resample_topography)
% 
% % Resample the model time step using the resample_topography variable
% dt = dt_topography;
% dt = dt * resample_topography;
% 
% 
% % Determine model length from first topography file
% filename = fullfile(path_model, 'topography.00000');
% topography = import_topography(filename);
% modellength = max(topography.x) / 1e3;
% 
% file = dir(fullfile(path_model, 'topography.*'));
% filenames = {file.name};
% num_files = numel(filenames);
% 
% % Create a new vector for interpolation with high resolution (10000 points)
% resolution_topo=10000;
% x_axis_interp = 0:360/resolution_topo:360;
% 
% 
% % Resample topo if needed for faster computing
% 
% % elevation_map = zeros(num_files, 10000 / (resample_topography * 1e3));
% i_start = 1;
% i_end = min(num_files, i_start + resample_topography - 1);
% 
% while i_start <= num_files
%     fprintf('Processing topography files %d of %d\n', i_start, num_files)
%     
%     for i = i_start:resample_topography:i_end
%         file_index = i;
% %         if i == 1
% %             file_index = 2;
% %         end
% 
%         % Import the data
%         topography = import_topography(fullfile(path_model, filenames{file_index}));
% 
%         % Sort the data in the desired order
%         x = topography.VarName1;
%         y = topography.x;
%         z = topography.y;
%         r = sqrt(x.^2 + y.^2);
%         %theta allows to know in which quadrant of the circle we are. 
%         theta = atan2(y, x); 
%         %now we can sort the data depending on theta index
%         [theta_sorted, theta_index] = sort(theta);
%         x_sorted = x(theta_index);
%         y_sorted = y(theta_index);
%         r_sorted = r(theta_index);
%         z_sorted = z(theta_index);
% 
%         % Convert theta to degrees and adjust the range to be 0 to 360 degrees
%         theta_degrees = theta_sorted * 180 / pi;
%         theta_degrees(theta_degrees<0) = theta_degrees(theta_degrees<0) + 360;
% 
%         % Find unique theta values and corresponding indices
%         [theta_unique, index_unique] = unique(theta_degrees);
% 
% 
%         % Interpolate the sorted z data to match the length of theta_degrees
%         z_interp = interp1(theta_unique, z_sorted(index_unique), x_axis_interp);
%         elevation_map(i, :) = z_interp ;
%     end
% 
%     % update index range for next iteration
%     i_start = i_end + 1;
%     i_end = min(num_files, i_start + resample_topography - 1);
% end
% 
% elevation = elevation_map(:, :);
% % time_elevation = 0:dt:dt * (size(elevation, 1) - 1);
% time_elevation=0:(dt*(num_files/resample_topography))/(num_files/resample_topography):dt*(num_files/resample_topography);
% end





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