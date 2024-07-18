function [time_var,average_params] = get_depth_averages(path_model, depths_for_average, averaged_parameters)
% Read data from text file
clear fileID;
filename = [path_model,'depth_average.txt'];
fileID = fopen(filename,'r');
header = fgetl(fileID); % Read the header line
header = strrep(header, '#', ''); % Remove '#' character from header
allowed_parameters = regexp(header, '\s+', 'split'); % Split header line using one or more spaces as delimiter
fclose(fileID);

allowed_parameters = allowed_parameters(~cellfun('isempty', allowed_parameters));


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

isDouble = cellfun(@(x) isnumeric(x) && all(cellfun(@isnumeric, num2cell(x))), dataArray);
dataArray = dataArray(isDouble);


parameter_indices = find(ismember(allowed_parameters, averaged_parameters));
parameter_indices = parameter_indices(1)-1; % Use only the first index

% Find index of depth parameter
depth_index = find(strcmp(allowed_parameters, 'depth'));
depth = dataArray{depth_index};

% Check if depths_for_average contains only depths that are in depth
available_depths = unique(depth);
if strcmp(depths_for_average, 'all')
    depths_for_average = available_depths';
    plot_all = 'true';
else
    plot_all = 'false';
    % Check if depths_for_average contains only depths that are in depth
    if ~isempty(setdiff(depths_for_average, available_depths))
        error(sprintf('Invalid depth(s) selected. The available depths are:\n%s', num2str(available_depths')));
    end
end

% Find the average values for each depth in depths_for_average
if ~strcmp(plot_all, 'true')
    time = dataArray{1}./1e6;
    average_params = cell(1, length(depths_for_average));
    for i=1:length(depths_for_average)
        d = depths_for_average(i);
        idx = depth == d;
        data_format= cell2mat(dataArray(parameter_indices));
        average_params{i} = data_format(idx);
        %     eval(averaged_parameters),'(idx)']);

    end
    time_var = time(idx);
end

if  strcmp(plot_all, 'true')
    % Create a table to hold all data
    average_params = table();
    for i = 1:length(allowed_parameters)
        average_params.(allowed_parameters{i}) = dataArray{i};
    end
    time_var = average_params.time./1e6;
end
end
