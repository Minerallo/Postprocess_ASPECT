function get_allstats(path_model,Output_numbers,cumultime,model_time_scaled,initial_state_time,path_shortening,data_number,shortening_stat,underthrusting_stat,color_time,colors,resampling_resolution,smoothing_stats_interval,plot_topography_time,dt,smoothing_interval,smooth_dip,resampling_topo,final_model_time,add_OP_vel)

cd(path_model);
% modellenght=2592;
% modelheight=900;
filename= 'statistics'; 
fid = fopen(filename,'r');
str = '#' ;
% Define how long is the head of the file and how many statistics
for c=1:999
        line=strtrim(fgets(fid)) ;
    if contains(line,str)
        stats_number = c ;
        %get each header take off semilicon to reveal the stat numbers and
        %title
        header(c)=string(line)
    else
        break; 
    end
end
fclose(fid) ;

%Get the data
for i =1:stats_number
    format_data= '%s'; 
    vector_format(i)=string(format_data);
end
final_format = join(vector_format); 
fid = fopen('statistics', 'rt');
Data = textscan(fid, final_format, 'headerLines', stats_number, 'CollectOutput', true);
fclose(fid);

%Convert the data to double
Data= str2double(Data{1,1});

%Plot all the data
Time = Data(:,2);

for i=1:stats_number
    figure(i);
    plot(Time,Data(:,i));xlabel('Time [yr]'); ylabel(header(i));title([header(i), 'versus Time']);  %Ma
end

% Plot a single graph defined
plot(Time,Data(:,data_number)+add_OP_vel);xlabel('Time [yr]'); ylabel(header(data_number));title([header(data_number), 'versus Time']);  %Ma

end

