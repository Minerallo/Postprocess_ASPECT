function [data_stats, stats_number, header, available_indices] = get_statistics(path_model, statistic_parameters)

    % Read the header of the statistics file
    filename = [path_model, 'statistics'];
    fid = fopen(filename, 'r');
    str = '#';
    stats_number = 0;
    header = string.empty;

    while ~feof(fid)
        line = strtrim(fgets(fid));
        if contains(line, str)
            % Found a new statistic header
            stats_number = stats_number + 1;
            header(stats_number) = line;
        end
    end

    fclose(fid);

    % Check which specified statistic parameters exist in the header
    available_parameters = cell(1, length(statistic_parameters));
    for i = 1:length(statistic_parameters)
        param_name = statistic_parameters{i};
        for j = 1:length(header)
            if contains(header(j), param_name)
                available_parameters{i} = header(j);
                break;
            end
        end
    end
    available_parameters(cellfun(@isempty, available_parameters)) = [];

    % Get the corresponding statistic numbers
    available_indices = [];
    for i = 1:length(available_parameters)
        available_indices(i) = str2double(regexp(extractBefore(available_parameters{i}, ':'), '\d+', 'match'));
    end

    if isempty(available_indices)
        % If no specified statistic parameters were found in the header
        error_msg = sprintf('The following specified statistic parameters were not found in the header:\n%s\nAvailable parameters: %s', strjoin(statistic_parameters, ', '), strjoin(header, ', '));
        error(error_msg);
    end

    %Get the data
    for i =1:stats_number
        format_data= '%s'; 
        vector_format(i)=string(format_data);
    end
    final_format = join(vector_format); 
    fid = fopen(filename, 'r');
    data_stats = textscan(fid, final_format, 'headerLines', stats_number, 'CollectOutput', true);
    fclose(fid);
    
    %Convert the data to double
    data_stats= str2double(data_stats{1,1});

end








% function [data,stats_number,header] = get_statistics(path, statistic_parameters)
% 
% filename = [path,'statistics'];
% fid = fopen(filename,'r');
% str = '#' ;
% % Define how long is the head of the file and how many statistics
% for c=1:999
%         line=strtrim(fgets(fid)) ;
%     if contains(line,str)
%         stats_number = c ;
%         %get each header take off semilicon to reveal the stat numbers and
%         %title
%         header(c)=string(line);
%     else
%         break; 
%     end
% end
% fclose(fid) ;
% 
% %Get the data
% for i =1:stats_number
%     format_data= '%s'; 
%     vector_format(i)=string(format_data);
% end
% final_format = join(vector_format); 
% fid = fopen('statistics', 'rt');
% data = textscan(fid, final_format, 'headerLines', stats_number, 'CollectOutput', true);
% fclose(fid);
% 
% %Convert the data to double
% data= str2double(data{1,1});
% 
% end

