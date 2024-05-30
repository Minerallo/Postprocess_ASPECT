function [init_step, max_step, end_step] = get_last_timestep(path_model,path_model_output,output_additional_maps_figures, visualize_model_at_specific_time,Interval_of_time_output_for_additional_postprocess)

num_files_surface=[];
num_files_lithosphere=[];
num_files_depths=[];

files_required_1 =[];
files_required_2 =[];
files_required_3 =[];

check_files_1 = [];
check_files_2 = [];
check_files_3 = [];
check_files_4 = [];
check_files_5 = [];
check_files_6 = [];
check_files_7 = [];

% Check if subduction_and_plume_number is true
% Surface files
file_surface = dir(fullfile(path_model, 'surface_*'));
filenames_surface = {file_surface.name};
num_files_surface = numel(filenames_surface);


file_depths = dir(fullfile(path_model, 'depths_*'));
filenames_depths = {file_depths.name};
num_files_depths = numel(filenames_depths);


file_lithosphere = dir(fullfile(path_model, 'lithosphere_*'));
filenames_lithosphere = {file_lithosphere.name};
num_files_lithosphere = numel(filenames_lithosphere);

max_num_files = max([num_files_lithosphere, num_files_surface, num_files_depths]);
max_step =max_num_files; 

for ite_check = 1:numel(output_additional_maps_figures)
check_files= dir(fullfile(path_model_output, strjoin(['*/' output_additional_maps_figures(ite_check) '_*'],"")));
    if isempty(check_files)
        num_check_files(ite_check) = 9999999;  % or any other default value you prefer
    else
        num_check_files(ite_check) = numel(check_files);
    end
num_check_files(ite_check) = numel(check_files);
end
min_check_files = min(num_check_files);

init_step=min_check_files;

if~isempty(visualize_model_at_specific_time)
    year_to_start = str2num(visualize_model_at_specific_time);
    end_step = round(year_to_start/Interval_of_time_output_for_additional_postprocess)+1;
else 
    end_step =[];
    year_to_start = [];
end
end

% % Assuming time and RMS_Surface_velocity are already defined
% % Assuming data_stats is a matrix containing time and velocity data
% 
% % Calculate RMS_Surface_velocity
% RMS_Surface_velocity = data_stats(:,39) .* data_stats(:,38);
% 
% % Plot the data with thicker line
% plot(time, RMS_Surface_velocity,'r', 'LineWidth', 2); % Set LineWidth to 2 for thicker line
% hold on; % Hold the plot to add legend and set axis limits
% 
% % Add legend
% % legend('RMS Surface Velocity');
% 
% % Set y-axis limits
% ylim([0 0.1]);
% 
% % Label the axes
% xlabel('Time');
% ylabel('RMS Surface Velocity');
% 
% % Add title
% title('RMS Surface Velocity over Time');

