function [time_elevation_layer, elevation_layer, x_axis_interp_layer, dip_layer] = get_topography_layer(path_model, dt_topography_layer, resample_topography_layer, model_length, model_height, calculate_topography_layer_dip, topography_smoothing_interval_for_dip_calculation,read_layer_elevation_from_bottom_to_top)

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

    % Check if the next iteration will go over num_files
    if (i > num_files)
        break;  % Exit the loop if the next iteration exceeds num_files
    end

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
    if strcmp(read_layer_elevation_from_bottom_to_top,'true')
        elevation_map(count_ite,:)=elevations2;
    else
        elevation_map(count_ite,:)=model_height-elevations2;
    end

%     if dip calculation is asked
    if strcmp(calculate_topography_layer_dip, 'true')
    
    % Check if the next iteration will go over num_files
    if i > num_files
        break;  % Exit the loop if the next iteration exceeds num_files
    end

%    find the indices of the places where there are no elevations2
%    (depracted from topolayer might not be relevant here)
            index_nan=find(isnan(elevations2));
%    The topography might be noisy because of the adaptative mesh refinement 
%    we can apply a slight smoothing to it
            smooth_topo=smooth(elevations2,topography_smoothing_interval_for_dip_calculation)';
%    Replace nan value by zero after smoothing (depracted from topolayer might not be relevant here)
            smooth_topo(index_nan)=0;
            dipsmooth=atand((diff(smooth_topo)./1e3));
            dip_map(count_ite,:)=dipsmooth;
%    For comparison with no smoothing uncomment the following lines
            elevations2(index_nan)=0;
            diptopo=atand((diff(elevations2)./1e3)); 
    end
end
elevation_layer = elevation_map(:, :);
time_elevation_layer = 0:dt:dt * (size(elevation_layer, 1) - 1);

if strcmp(calculate_topography_layer_dip, 'true')
    dip_layer = dip_map(:, :);
%    The figure will plot the last profil for comparison with and without smoothing
    figure();
    plot(1:numel(diptopo), diptopo, 'r-', 'LineWidth', 1.4);
    hold on;
    plot(1:numel(dipsmooth), dipsmooth, 'b-', 'LineWidth', 1.4);
    legend('Dip Not Smoothed','Dip Smoothed');
    xlabel('Model length');
    ylabel('Dip °');title('Dip Layer Comparison');legend('Location', 'Best'); 
else
    dip_layer = 0;
end

end