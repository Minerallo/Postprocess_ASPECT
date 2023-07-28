function [time_elevation_layer, elevation_layer, x_axis_interp_layer] = get_topography_layer(path_model, dt_topography_layer, resample_topography_layer, model_length, model_height)

% Resample the model time step using the resample_topography_layer variable
dt = dt_topography_layer;
dt = dt * resample_topography_layer;


% Determine model length from first topography file

file = dir(fullfile(path_model, 'topolayer.*'));
filenames = {file.name};
num_files = numel(filenames);

% for i = 1:resample_topography_layer:num_files 
% filename_local = fullfile(path_model, filenames{i});
% topography_local = import_topography_layer(filename_local);
% % Give the highest x value of the file
% modellength_local = max(topography.VarName1) / 1e3;
% % Give the highest x value in any of the file
% if (modellength_local>modellength)
% modellength =modellength_local ; 
% end
% end

% Create a new vector for interpolation with high resolution (10000 points)
x_axis_interp_layer = 0:1e3:model_length;

% Resample topo if needed for faster computing
count_ite=0;
% elevation_map = zeros(num_files, 10000 / (resample_topography_layer * 1e3));
for i = 1:resample_topography_layer:num_files
    count_ite=count_ite+1;
    fprintf('Processing topography layer file %d of %d\n', i, num_files)
    file_index = i;
%     if i == 1
%         file_index = 2;
%     end
    
    % Import the data
    topolayer = import_topography_layer(fullfile(path_model, filenames{file_index}));

    %% Convert to output type
    topolayer = table2array(topolayer);
    
    maxdepth=max(max(topolayer(:,2)));
    mindepth=min(min(topolayer(:,2)));
    [val_sorted index_sorted]= sortrows(topolayer,[1 2]);
    [Xunique_sorted, Xunique_index] = unique(val_sorted(:,1),'last');
    elevations1=val_sorted(Xunique_index,2);
    elevations2=interp1(Xunique_sorted,elevations1,x_axis_interp_layer);
    elevation_map(count_ite,:)=model_height-elevations2;
    
end
elevation_layer = elevation_map(:, :);
time_elevation_layer = 0:dt:dt * (size(elevation_layer, 1) - 1);
end