function [time_elevation_topography, elevation_topography, x_axis_interp, dip_topography] = get_topography_box(path_model, dt_topography, resample_topography, model_length, model_height, calculate_topography_dip, topography_smoothing_interval_for_dip_calculation)


% Resample the model time step using the resample_topography variable
dt = dt_topography;
dt = dt * resample_topography;


% Determine model length from first topography file

file = dir(fullfile(path_model, 'topography.*'));
filenames = {file.name};
num_files = numel(filenames);

% Create a new vector for interpolation with high resolution (10000 points)
x_axis_interp = 0:1e3:model_length;

% Resample topo if needed for faster computing
count_ite=0;
% elevation_map = zeros(num_files, 10000 / (resample_topography * 1e3));
for i = 1:resample_topography:num_files
    count_ite=count_ite+1;

    % Check if the next iteration will go over num_files
    if (i > num_files)
        break;  % Exit the loop if the next iteration exceeds num_files
    end

    fprintf('Processing topography box file %d of %d\n', i, num_files)
    file_index = i;
%     if i == 1
%         file_index = 2;
%     end
    
    % Import the data
    topobox = import_topography_box(fullfile(path_model, filenames{file_index}));

    %% Convert to output type
%     topobox = table2array(topobox);
    [val_sorted index_sorted]= sortrows(topobox,[1 2]);
    [Xunique_sorted, Xunique_index] = unique(val_sorted.x,'last');

%     elevations1=val_sorted(Xunique_index,2);
    elevations2=interp1(Xunique_sorted,topobox.h(Xunique_index),x_axis_interp');
%         elevation_map(count_ite,:)=elevations2;

        elevation_map(count_ite,:)=elevations2;


%     if dip calculation is asked
    if strcmp(calculate_topography_dip, 'true')
    
    % Check if the next iteration will go over num_files
    if i > num_files
        break;  % Exit the loop if the next iteration exceeds num_files
    end

%    find the indices of the places where there are no elevations2
%    (depracted from topobox might not be relevant here)
            index_nan=find(isnan(elevations2));
%    The topography might be noisy because of the adaptative mesh refinement 
%    we can apply a slight smoothing to it
            smooth_topo=smooth(elevations2,topography_smoothing_interval_for_dip_calculation)';
%    Replace nan value by zero after smoothing (depracted from topobox might not be relevant here)
            smooth_topo(index_nan)=0;
            dipsmooth=atand((diff(smooth_topo)./1e3));
            dip_map(count_ite,:)=dipsmooth;
%    For comparison with no smoothing uncomment the following lines
            elevations2(index_nan)=0;
            diptopo=atand((diff(elevations2)./1e3)); 
    end
end
elevation_topography = elevation_map(:, :);
time_elevation_topography = 0:dt:dt * (size(elevation_topography, 1) - 1);

if strcmp(calculate_topography_dip, 'true')
    dip_topography = dip_map(:, :);

else
    dip_topography = 0;
end

end