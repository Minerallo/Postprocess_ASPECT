function [time_avg,param_avg] = get_depth_averages(path_model, depths_for_average, averaged_parameters)
% Read data from text file
clear fileID;
filename = [path_model,'depth_average.txt'];
fileID = fopen(filename,'r');
header = fgetl(fileID); % Read the header line
header = strrep(header, '#', ''); % Remove '#' character from header
allowed_parameters = regexp(header, '\s+', 'split'); % Split header line using one or more spaces as delimiter
fclose(fileID);

% Check if averaged_parameter is allowed
if ~ismember(averaged_parameters, allowed_parameters)
    error(sprintf('Invalid parameter selected. The allowed parameters are:\n%s', strjoin(allowed_parameters, ', ')));
end

% Read data from text file
delimiter = ' ';
formatSpec = ['%f' repmat('%f',1,length(allowed_parameters)-1) '%[^\n\r]'];
fileID = fopen(filename,'r');
dataArray = textscan(fileID, formatSpec, 'Delimiter', delimiter, 'MultipleDelimsAsOne', true, 'HeaderLines' ,1, 'ReturnOnError', false);
fclose(fileID);


parameter_indices = find(ismember(allowed_parameters, averaged_parameters));
parameter_indices = parameter_indices(1)-1; % Use only the first index


% Find index of depth parameter
depth_index = find(strcmp(allowed_parameters, 'depth'));
depth_index = depth_index -1;
depth = dataArray{depth_index}; 

% Check if depths_for_average contains only depths that are in depth
available_depths = unique(depth);
if ~isempty(setdiff(depths_for_average, available_depths))
    error(sprintf('Invalid depth(s) selected. The available depths are:\n%s', num2str(available_depths')));
end

% Find the average values for each depth in depths_for_average
time = dataArray{1}./1e6;
param_avg = cell(1, length(depths_for_average));
for i=1:length(depths_for_average)
    d = depths_for_average(i);
    idx = depth == d;
    data_format= cell2mat(dataArray(parameter_indices));
    param_avg{i} = data_format(idx);
%     eval(averaged_parameters),'(idx)']);

end
time_avg = time(idx);

end
