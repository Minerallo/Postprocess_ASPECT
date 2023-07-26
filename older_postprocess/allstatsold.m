function allstats(path_model,Output_numbers,cumultime,model_time_scaled,initial_state_time,path_shortening,data_number,shortening_stat,underthrusting_stat,color_time,colors,resampling_resolution,smoothing_stats_interval,plot_topography_time,dt,smoothing_interval,smooth_dip,resampling_topo,final_model_time,add_OP_vel)
cd(path_model);
% data_number=77; 
modellenght=2592;
modelheight=900;
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

% for i=1:stats_number
%     figure(i);
%     plot(Time,Data(:,i));xlabel('Time [yr]'); ylabel(header(i));title([header(i), 'versus Time']);  %Ma
% end

% Plot a single graph defined
% plot(Time,Data(:,data_number)+add_OP_vel);xlabel('Time [yr]'); ylabel(header(data_number));title([header(data_number), 'versus Time']);  %Ma



%% Velocity module
if any(Output_numbers==[8,5])
%%To compare to Sdrolias 2006 
%Plot convergence rate of Sdrolias for comparison
Sdrolias=[20, 14.524625267665952
15.005096839959226, 14.524625267665952
14.984709480122325, 12.537473233404711
10.010193679918451, 12.50321199143469
9.96941896024465, 8.289079229122056
5.03567787971458, 8.289079229122056
4.974515800203875, 7.809421841541756
0, 7.809421841541756];
abs_downgoing_plate_velocity = [2.4283018169707304, 51.194661337872645
7.420345631499089, 83.1335035586969
12.427191651556072, 94.41093389581712
17.43181734078378, 108.78757601549296
22.25659623283869, 124.20037252217244
27.45661103504422, 115.84637777695542
32.32431632313214, 71.34107982089327
37.522850904784804, 65.05322626404666
42.32246604744108, 115.59042297302295
47.16944824778892, 100.01110165414647
52.385005365799515, 69.96262443104018
57.6101839174037, 26.484229483526377];


% line_width=1.5;
% figure(1);clf;close all;
% set(gcf,'color','w');  %Ma
% plot((flip((abs_downgoing_plate_velocity(:,1)-32.5)).*1e6-30e6)./1e6,abs_downgoing_plate_velocity(:,2)./1e3,'ro--','LineWidth',line_width);hold on;xlabel('Time [Ma]'); ylabel("Oceanic Plate Velocity [cm/yr]");title("Oceanic plate velocity");
% plot(Time./1e6-30,Data(:,data_number)+add_OP_vel,'b-','LineWidth',line_width);
% legend("Paleomagnetic (Sdrolias et al., 2006)","Model")

time_rescaled =model_time_scaled./1e6; %rescaled to geological time
line_width=1.5;
figure(5);clf;

time_model_for_plot=Time./1e6;
for c=1:numel(color_time)
[value_residual,index_color{c}] = min(abs(time_model_for_plot.*1e6 - color_time(c)));
end


startStopIdx = [1 index_color{1};    
     index_color{1} index_color{2};
    index_color{2} index_color{3}; 
    index_color{3} index_color{4};
    index_color{4} numel(time_model_for_plot)]; % [start, stop]



timing_vel=time_model_for_plot;
%oroginal poster color
% colors= [233 235 246; 245 238 245 ; 253 234 236;254 253 237]; 
%new color underthruting to flatslab time respectively
startStopX = timing_vel(flip(startStopIdx))-model_time_scaled/1e6; 
width = startStopX(:,2)-startStopX(:,1);  

yli(1)=0.001;
yli(2)=0.15-0.01;


arrayfun(@(i)rectangle('Position', [startStopX(i,1),yli(1),width(i),range(yli)], ...
    'EdgeColor', colors(i,:), 'FaceColor', colors(i,:)), 1:size(startStopX,1)); hold on;

set(gcf,'color','w');  %Ma
plot(flip(abs_downgoing_plate_velocity(:,1))-abs_downgoing_plate_velocity(end,1),abs_downgoing_plate_velocity(:,2)./1e3,'ro--','LineWidth',line_width);hold on;xlabel('Time [Ma]'); ylabel("Oceanic Plate Velocity [m/yr]");title("Oceanic plate velocity");
plot(Time./1e6-time_rescaled,Data(:,data_number)+add_OP_vel,'b-','LineWidth',line_width);
legend("Paleomagnetic (Sdrolias et al., 2006)","Model",'Location','northwest');


figure(8);clf;
vel = Data(:,data_number); 
acceleration=diff(vel)./diff(Time);
acceleration=smooth(acceleration,100);
% yli(1)=0.001;
% yli(2)=0.15-0.01;
% 
% 
% arrayfun(@(i)rectangle('Position', [startStopX(i,1),yli(1),width(i),range(yli)], ...
%     'EdgeColor', colors(i,:), 'FaceColor', colors(i,:)), 1:size(startStopX,1)); hold on;

set(gcf,'color','w');  %Ma
plot(Time(2:end)./1e6-time_rescaled,acceleration,'b-','LineWidth',line_width);
legend("Modelled slab acceleration",'Location','northwest');ylim([-1e-8 1e-8]); 

end

%% Shortening module
if any(Output_numbers==[8,6,7])
    cd(path_shortening)
%Plot all the data
load('shortening_rate_oncken.mat')
load('oncken2012.mat');

Time = Data(:,2);
IRS=6; %initial refinement steps number
% initial_state_time = 12e6; %shortening difference will be calculated from that time , write 1 to start from 0Ma
% shortening_stat=84;
% underthrusting_stat=85;
initial_state_null=true;

%Find the reference index from which the shortenning will be calculated
ind_ref = find(min(abs(Time - initial_state_time)) == abs(Time - initial_state_time),1,'last');

%Calculate the shortening from reference time given
shortening_vec=Data(ind_ref,shortening_stat)-Data(:,shortening_stat);
shortening_foreland_vec=Data(:,underthrusting_stat)-Data(ind_ref,underthrusting_stat);

%%Avoiding IRS steps for the interpolation
Time = Time(IRS:end);
shortening_vec=shortening_vec(IRS:end);
shortening_foreland_vec=shortening_foreland_vec(IRS:end);

%%Something is off yith the following but I let it there for later
%%if the initial state time is is not 0Ma (=1) then cut the values of the initial state
%Otherwise comment these 3 lines or give initial state time a value of 1
% Time = Time(ind_ref:end);
% shortening_vec=shortening_vec(ind_ref:end);
% shortening_foreland_vec=shortening_foreland_vec(ind_ref:end);

%%Delete the nan values at the end 

%Resample to lower resolution so it will help for the smoothing later
time_vec=0:resampling_resolution:model_time_scaled;
shortening_vec=interp1(Time,shortening_vec,time_vec);
shortening_foreland_vec=interp1(Time,shortening_foreland_vec,time_vec);

%%Get rid of the nan values so they dont interfere with the smoothing later
idx = arrayfun(@(v)isnumeric(v)&&any(isnan(v)),shortening_vec);
shortening_vec= shortening_vec(1:numel(shortening_vec)-nnz(idx));
shortening_foreland_vec= shortening_foreland_vec(1:numel(shortening_foreland_vec)-nnz(idx));
time_vec= time_vec(1:numel(time_vec)-nnz(idx));

% figure(7), plot(orogen_shortening_rate);hold on; plot(orogen_shortening_rate_smooth)


%Calculate the shortening rate
time_interval=diff(time_vec); 
orogen_shortening_rate=diff(shortening_vec.*1e3)./time_interval;
underthrusting_rate=diff(shortening_foreland_vec.*1e3)./time_interval;



if (initial_state_null)

  [minValue,closestIndex] = min(abs(time_vec - initial_state_time));
  orogen_shortening_rate(1:closestIndex-1)=0;
  underthrusting_rate(1:closestIndex-1)=0;
end

% time_to_plot =(model_time_scaled/1e6)-time_to_plot1;
% time_to_plot =(model_time_scaled/1e6)-time_to_plot1;
time_vec_to_plot=model_time_scaled/1e6-time_vec./1e6;

total_shortening = shortening_vec +shortening_foreland_vec;
total_shortening_rate = orogen_shortening_rate  +underthrusting_rate;



% oncken_time=0:55; 
% model_time = 46_

%Time scale for the shortening rate
cum_time_to_plot1=cumsum(time_interval)./1e6;
cum_time_to_plot =model_time_scaled/1e6-cum_time_to_plot1;

%Smoothing of the shortening rate
total_shortening_rate_smooth=smooth(total_shortening_rate,smoothing_stats_interval);
orogen_shortening_rate_smooth=smooth(orogen_shortening_rate,smoothing_stats_interval);
underthrusting_shortening_rate_smooth=smooth(underthrusting_rate,smoothing_stats_interval);
end

% figure(7);clf;
% set(gca, 'XDir','reverse');hold on;
% % plot(cum_time_to_plot,total_shortening_rate,'LineWidth',2); 
% plot(cum_time_to_plot,total_shortening_rate_smooth,'LineWidth',2);
% % plot(time_to_plot,orogen_shortening_rate,'LineWidth',2);
% % plot(time_to_plot,underthrusting_shortening_rate,'LineWidth',2);
% plot(cum_time_to_plot,orogen_shortening_rate_smooth,'LineWidth',2);
% plot(cum_time_to_plot,underthrusting_shortening_rate_smooth,'LineWidth',2);
% plot(0:55,cumulative_shortening_Ma,'r-*','LineWidth',2);
% plot(oncken2012(:,1),oncken2012(:,2),'k-*','LineWidth',2)
% title ('Total shortening rate evolution of Southern Central Andes (21S)');xlabel('Geological time [Ma]');ylabel('shortening rate [mm/a]');set(gcf,'color','w');
% legend('Total shortening rate s5Ma','Orogenic shortening rate s5Ma', 'Underthrusting shortening rate s5Ma','Oncken 2006','Oncken 2012','Location','northwest');
% % legend('total shortening rate','total shortening rate s5Ma','orogen shortening rate s5Ma', 'underthrust shortening rate s5Ma','Oncken 2006','Oncken 2012','Location','northwest');
% ylim([-10 20]);xlim([0 60]);

if any(Output_numbers==[8,6])
maxtime=50;
fig=figure(7);clf;
% setup bottom axis
ax = axes();
hold(ax);
ax.YAxis.Scale = 'linear';
xlabel(ax, 'Geological time [Ma]', 'Interpreter', 'latex', 'FontSize', 14,'fontweight','bold');
ylabel(ax, 'shortening rate [mm/a]', 'Interpreter', 'latex', 'FontSize', 14,'fontweight','bold');
% setup top axis
ax_top = axes('Parent', fig); % axis to appear at top
ax_top.Position = ax.Position;
ax_top.YAxis.Visible = 'off';
ax_top.XAxisLocation = 'top';
% ax_top.XDir = 'reverse';
ax_top.Color = 'none';
xlabel(ax_top, 'Model time [Ma]', 'Interpreter', 'latex', 'FontSize', 14,'fontweight','bold');
ax.XLim = [-maxtime 0];ax.YLim = [-10 20];
ax_top.XLim = [0 maxtime];
y_ticks = [model_time_scaled./1e6-(abs(maxtime-50)) model_time_scaled./1e6-(abs(maxtime-45)) model_time_scaled./1e6-(abs(maxtime-40)) model_time_scaled./1e6-(abs(maxtime-35)) model_time_scaled./1e6-(abs(maxtime-30)) model_time_scaled./1e6-(abs(maxtime-25)) model_time_scaled./1e6-(abs(maxtime-20)) model_time_scaled./1e6-(abs(maxtime-15)) model_time_scaled./1e6-(abs(maxtime-10)) 5-(maxtime-model_time_scaled./1e6) -6];

% % find phases of deformation
% %if it was not done before we need to correct the initialization time that
% %send back some high shortening rate value in some case
%   [minValue,closestIndex] = min(abs(flip(cum_time_to_plot.*1e6) - initial_state_time));
%   total_shortening_rate_smooth(1:closestIndex-1)=0;
%   orogen_shortening_rate_smooth(1:closestIndex-1)=0;
%   underthrusting_shortening_rate_smooth(1:closestIndex-1)=0;
%   %Give the index close to 20 Ma
%   [value_residual,index_close_to] = min(abs(flip(cum_time_to_plot.*1e6) - (initial_state_time+4e6)));
%   %I cut the vector at 20 Ma so I dont catch the steepening phase
%   buckling_vector_considered = orogen_shortening_rate_smooth(index_close_to:end);
% % Track underthrusting index  
% [ind_underthrusting]=find(underthrusting_shortening_rate_smooth>=1.5);
% diff_buckling_vector= diff(buckling_vector_considered); 
% [ind_buckling]=find(diff_buckling_vector>=1,1);
% %Take 2 index before maximum of buckling shortening rate
% ind_buckling=index_close_to+ind_buckling; 
% 
% %adjust to time from window
% 
% max_x_plot = ax_top.XLim(2)*1e6
% difference_time_model=max_x_plot-model_time_scaled
% init_time_value=-(max_x_plot-difference_time_model)
% 
% if(isempty(ind_underthrusting))
%     ind_min_underthrusting = numel(cum_time_to_plot);
%      ind_max_underthrusting = numel(cum_time_to_plot);
% else
%         ind_min_underthrusting = min(ind_underthrusting);
%      ind_max_underthrusting = max(ind_underthrusting);
% end
% 
% if(isempty(ind_buckling))
%     ind_buckling = ind_min_underthrusting;
% end
% % Specify start, stop index values (positive integers)
% startStopIdx = [1 closestIndex;       
%     closestIndex ind_buckling;
%     ind_buckling ind_min_underthrusting; 
%     ind_min_underthrusting ind_max_underthrusting]; % [start, stop]

for c=1:numel(color_time)
[value_residual,index_color{c}] = min(abs(flip(cum_time_to_plot.*1e6) - color_time(c)));
end

startStopIdx = [1 index_color{1};       
    index_color{1} index_color{2};
    index_color{2} index_color{3}; 
    index_color{3} index_color{4};
    index_color{4} numel(cum_time_to_plot)]; % [start, stop]
max_x_plot = ax_top.XLim(2)*1e6;
difference_time_model=max_x_plot-model_time_scaled;


% Convert the index values to rectangle coordinates
ylimi=ax.YLim;
ylimi(1)=ylimi(1);
ylimi(2)=ylimi(2);
ylim([ylimi(1)-0.05 ylimi(2)+0.1]); 
yl = ylimi;
timing=flip(cum_time_to_plot);
%oroginal poster color
%  colors= [233 235 246; 245 238 245 ; 253 234 236;254 253 237]; 
%new color underthruting to flatslab time respectively

startStopX = timing(flip(startStopIdx))+difference_time_model/1e6+0.5; 
width = startStopX(:,2)-startStopX(:,1);  


hold(ax_top,'on');

arrayfun(@(i)rectangle('Parent',ax_top,'Position', [startStopX(i,1),yl(1),width(i),range(yl)], ...
     'EdgeColor', colors(i,:),'FaceColor', colors(i,:)), 1:size(startStopX,1));hold on;%,''EdgeColor', 'black', 'LineWidth',0.5,LineStyle,':''


% second_xaxes = model_time_scaled./1e6-y_ticks;
% ax_top.XTick = fliplr(second_xaxes);
ax_top.XTickLabel = compose('%1.0f',fliplr(y_ticks));
% plot(ax,cum_time_to_plot,total_shortening_rate,'LineWidth',2); 
 plot(ax_top,ax_top.XLim(2)-cum_time_to_plot,total_shortening_rate_smooth,'LineWidth',2);
% plot(ax,time_to_plot,orogen_shortening_rate,'LineWidth',2);
% plot(ax,time_to_plot,underthrusting_shortening_rate,'LineWidth',2);
plot(ax_top,ax_top.XLim(2)-cum_time_to_plot,orogen_shortening_rate_smooth,'LineWidth',2);
plot(ax_top,ax_top.XLim(2)-cum_time_to_plot,underthrusting_shortening_rate_smooth,'LineWidth',2);
plot(ax_top,ax_top.XLim(2)-(0:55),cumulative_shortening_Ma,'r--','LineWidth',2);
plot(ax_top,ax_top.XLim(2)-abs(oncken2012(end,1)-oncken2012(:,1)),oncken2012(:,2),'k-','LineWidth',2);set(gcf,'color','w');
title ('Comparison between the modelled and the observed shortening rate at 21°S');set(gcf,'color','w');
legend(ax_top,'Total shortening rate (model)','Orogenic shortening rate (model)', 'Underthrusting shortening rate (model)','Geological data Oncken 2006','Geological data Oncken 2012','Location','northwest','FontSize',10);
wind=get(gca,'position'); set(gca,'position',wind);set(gca, 'FontName', 'Arial');
% fig.Position = [0 0 5120 2880];
%saveas(gcf,[path_model '.png']); 
hold off ;
end

if any(Output_numbers==[8,6])
figdos=figure(6);clf;
ax6= axes('Parent', figdos);
  [minValue,closestIndexfull] = min(abs(flip(time_vec_to_plot.*1e6) - initial_state_time));
  total_shortening(1:closestIndexfull-1)=0;
  shortening_vec(1:closestIndexfull-1)=0;
  shortening_foreland_vec(1:closestIndexfull-1)=0;
  
%   % Specify start, stop index values (positive integers)
% startStopIdx = [1 closestIndexfull;       
%     closestIndexfull ind_buckling;
%     ind_buckling ind_min_underthrusting; 
%     ind_min_underthrusting ind_max_underthrusting]; % [start, stop]
% % Convert the index values to rectangle coordinates
% 

for c=1:numel(color_time)
[value_residual,index_color{c}] = min(abs(flip(time_vec_to_plot.*1e6) - color_time(c)));
end


startStopIdx = [1 index_color{1};    
     index_color{1} index_color{2};
    index_color{2} index_color{3}; 
    index_color{3} index_color{4};
    index_color{4} numel(time_vec_to_plot)]; % [start, stop]



timing_short=flip(time_vec_to_plot);
%oroginal poster color
% colors= [233 235 246; 245 238 245 ; 253 234 236;254 253 237]; 
%new color underthruting to flatslab time respectively
startStopX = timing_short(flip(startStopIdx)); 
startStopX(end,1)=1;
width = startStopX(:,2)-startStopX(:,1);  

yli(1)=-15+1;
yli(2)=350-1;


arrayfun(@(i)rectangle('Position', [startStopX(i,1),yli(1),width(i),range(yli)], ...
    'EdgeColor', colors(i,:), 'FaceColor', colors(i,:)), 1:size(startStopX,1)); hold on;

plot(max(time_vec_to_plot)-time_vec_to_plot, total_shortening./1e3,'LineWidth',2)
plot(max(time_vec_to_plot)- time_vec_to_plot, shortening_vec./1e3,'LineWidth',2);%set(gca, 'XDir','reverse');
plot(max(time_vec_to_plot)-time_vec_to_plot, shortening_foreland_vec./1e3,'LineWidth',2);%set(gca, 'XDir','reverse');
title ('Modelled shortening of Central Andes 21°S');xlabel('Model time [My]');ylabel('shortening [km]');set(gcf,'color','w');
legend('Shortening total','Shortening orogen', 'shortening underthrusting','Location','northwest' );ylim([-15 350]); xlim([0 model_time_scaled/1e6]);

end

if any(Output_numbers==[7,8])
shortening_at_t=cumsum(total_shortening_rate_smooth);
t_want=cumultime./1e6; % give model time
track_t=abs(t_want-cum_time_to_plot);
[indmin_t]=find(track_t==min(abs(t_want-cum_time_to_plot)));
%get cumulative shortening
shortening_at_t(indmin_t) 
end


dt=dt*resampling_topo;
%topography plot using ASPECT postprocess
cd(path_model);
file = dir ('topography.*');
filenames = {file.name};
A=numel(filenames);

%A variable in case of quick resampling
z=0;
if any(Output_numbers==[1,2,3,4,8])
    for i=0:resampling_topo:A
        file_index=i;
        if (i==0)
            file_index=1;
        end
        z=z+1;
        %resample topo by the number given for faster computing
        pathelev = sprintf('%s',pwd,'/',filenames{file_index});
        fprintf('%i/%i\n', file_index,A)
        %% Setup the Import Options and import the data
        opts = delimitedTextImportOptions("NumVariables", 4);
        
        % Specify range and delimiter
        opts.DataLines = [2, Inf];
        opts.Delimiter = " ";
        
        % Specify column names and types
        opts.VariableNames = ["VarName1", "x", "y", "Var4"];
        opts.SelectedVariableNames = ["VarName1", "x", "y"];
        opts.VariableTypes = ["double", "double", "double", "string"];
        
        % Specify file level properties
        opts.ExtraColumnsRule = "ignore";
        opts.EmptyLineRule = "read";
        opts.ConsecutiveDelimitersRule = "join";
        opts.LeadingDelimitersRule = "ignore";
        
        % Specify variable properties
        opts = setvaropts(opts, "Var4", "WhitespaceRule", "preserve");
        opts = setvaropts(opts, "Var4", "EmptyFieldRule", "auto");
        
        % Import the data
        topography = readtable(pathelev, opts);
        
        %% Convert to output type
        topography = table2array(topography);
        
        maxdepth=max(max(topography(:,3)));
        
        [Xsorted,ind] = sort(topography(:,1));
        elevations1=topography(ind,3);
        [Xunique, index] = unique(Xsorted);
        elevations2=interp1(Xunique,elevations1(index),0:1e3:modellenght*1e3);
        Elevationmap(z,:)=elevations2;
    end
end
line_width =2;
time=0:final_model_time/(A/resampling_topo):final_model_time
%     time=0:dt:(numel(1:resampling_topo:A)-1)*dt;
bathymetry=3e3;

if any(Output_numbers==[1,3,8])
    figure(1);clf;
    surf(1:size(Elevationmap,2),time,Elevationmap./1e3-bathymetry./1e3);shading interp;c=colorbar;demcmap('inc',[5 -8],0.1);ylabel('Time[Ma]'),xlabel('Length[km]');zlabel('Elevation[km]');ylim([0 model_time_scaled]);xlim([0 modellenght]);set(gcf,'color','w');
    %  surf(1:size(Elevationmap,2),time,Elevationmap./1e3-bathymetry./1e3);shading interp;c=colorbar;demcmap('inc',[(max(max(Elevationmap))-Ridge)/1e3 (min(min(Elevationmap))-Ridge)/1e3],0.1);ylabel('Time[Ma]'),xlabel('Length[km]');zlabel('Elevation[km]');
    c.Label.String= "Elevations [km]";set(gcf,'color','w');view(2);%set(gca, 'color', 'none');grid off;set(gca,'XColor', 'none','YColor','none','ZColor','none');
    % % FaceLighting = 'gour
end



if any(Output_numbers==[1,8])
    %comparison with real topography
    
    %Import topography profil(extracted from geomappapp)
    opts = delimitedTextImportOptions("NumVariables", 4);
    opts.DataLines = [3, Inf];
    opts.Delimiter = "\t";
    opts.VariableNames = ["VarName1", "VarName2", "VarName3", "GMRTGridVersion37"];
    opts.VariableTypes = ["double", "double", "double", "double"];
    opts.ExtraColumnsRule = "ignore";
    opts.EmptyLineRule = "read";
    topogreatcircle = readtable("/home/ponsm/Nextcloud/phd-central-andes-shortening/cookbookmika/Plot_topography/topo21_greatcircle.txt", opts);
    topogreatcircle = table2array(topogreatcircle);
    clear opts
    
    %To georeference then lets use the minimum elevation, the trench.
    indexmin_model=find(Elevationmap(end,:)==min(Elevationmap(end,:)));
    indexmin_topo =find(topogreatcircle(:,end)==min(topogreatcircle(:,end)));
    Vec_lenght= 1:size(Elevationmap,2);
    
    % Correction from water
    
    topo=topogreatcircle(:,end);
    % [xx,yy]=find(topo==min(topo));
    xx=find(topo>0,1);
    sea_level_correction=topo(1:xx)+(1030/3300)*abs(topo(1:xx));
    topocorrected(1:xx) = sea_level_correction(1:xx);
    topocorrected(xx+1:size(topo,1)) = topo(xx+1:end);
    
    figure(2); clf;
    plot(topogreatcircle(:,3)-topogreatcircle(indexmin_topo,3),topogreatcircle(:,end)); hold on ;
    plot(topogreatcircle(:,3)-topogreatcircle(indexmin_topo,3),topocorrected);
    plot(Vec_lenght-Vec_lenght(indexmin_model(1)),Elevationmap(end,:)-bathymetry);
    legend('Real topography','Topography corrected from water','Model');xlabel('Distance from the trench [km]');
    ylabel('Altitude [km]'); title('Topographic comparison'); set(gcf,'color','w');%set(gca, 'color', 'none');
end


%% %%%%%%%%%%
if any(Output_numbers==[1,3,8])
    %Plot every 10 Ma
    dt_elevation=10e6/dt;
    
    figure(3); clf;
    % plot(topogreatcircle(:,3)-topogreatcircle(indexmin_topo,3),topogreatcircle(:,end));
    plot(topogreatcircle(:,3)-topogreatcircle(indexmin_topo,3),topocorrected/1e3,'LineWidth',line_width);hold on ;
    for i=1:fix(size(Elevationmap,1)/dt_elevation)
        plot(Vec_lenght-Vec_lenght(indexmin_model(1)),(Elevationmap(i*dt_elevation,:)-bathymetry)/1e3,'LineWidth',line_width);
        
    end
    xlabel('Distance from the trench [km]');
    ylabel('Altitude [km]'); title('Evolution of surface topography'); set(gcf,'color','w');%set(gca, 'color', 'none');
    legend ('Topography corrected from water','10 Ma','20 Ma','30 Ma','40 Ma','50 Ma');
end

% figure(3);clf;
% plot(topogreatcircle(:,3)-topogreatcircle(indexmin_topo,3),topocorrected/1e3,'LineWidth',line_width);hold on ;
% plot(Vec_lenght-Vec_lenght(indexmin_model(1)),(Elevationmap(end,:)-bathymetry)/1e3,'LineWidth',line_width);
% legend('Real topography corrected from water','Model');set(gcf,'color','w');
% xlabel('Distance from the trench [km]');

if any(Output_numbers==[1,4,8])
    %% trench migration rate
    time_to_plot=model_time_scaled./1e6-time./1e6;
    load('/home/ponsm/Nextcloud/phd-central-andes-shortening/cookbookmika/Plot_topography/trench_rollback.mat')
    for i=2:size(Elevationmap,1)
        indexmin_topo_new(i)=find(Elevationmap(i,:)==min(Elevationmap(i,:)),1);
    end
    % for beginning of the model topography isn't right replace retreat by 0
    trench_rate=diff(indexmin_topo_new).*1e5./dt*-1; %-1 by convention to have westward motion positive
    idx_gap=find(abs(trench_rate)>50);
    for a=1:numel(idx_gap)
        trench_rate(idx_gap(a))=0;
    end
    trench_rate=[0 trench_rate] ;% add a 0 in front so it fit with the length of time
    trench_rate_smooth=smooth(trench_rate,smoothing_interval);
    
    %Do not comptabilize before initial time
    [minValue,closestIndexfull] = min(abs(time - initial_state_time));
    trench_rate(1:closestIndexfull-1)=0;
    trench_rate_smooth(1:closestIndexfull-1)=0;
    
    figure(4); clf;
    time_for_colors=flip(time_to_plot);
    for c=1:numel(color_time)
        [value_residual,index_color{c}] = min(abs(time_for_colors*1e6 - color_time(c)));
    end
    
    
    startStopIdx = [1 index_color{1};
        index_color{1} index_color{2};
        index_color{2} index_color{3};
        index_color{3} index_color{4};
        index_color{4} numel(time_for_colors)]; % [start, stop]
        startStopX = time_for_colors(flip(startStopIdx));
    startStopX(end,1)=0.2;
    width = startStopX(:,2)-startStopX(:,1);
    
    
    yli(1)=-2+0.5;
    yli(2)=4-0.5;
    
    
    arrayfun(@(i)rectangle('Position', [startStopX(i,1),yli(1),width(i),range(yli)], ...
        'EdgeColor', colors(i,:), 'FaceColor', colors(i,:)), 1:size(startStopX,1)); hold on;
    
    plot(max(time_to_plot)-time_to_plot,trench_rate,'LineWidth',2);hold on;
    plot(max(time_to_plot)-time_to_plot,trench_rate_smooth,'LineWidth',2);
    plot(model_time_scaled/1e6-trench_rollback(:,1),trench_rollback(:,2),'LineWidth',2)
    % plot(trench_rollback(:,1),trench_rollback(:,2)*10,'LineWidth',2)
    legend('trench velocity','trench velocity smoothed 5Ma','Reconstruction Oncken2006');
    ylabel('trench migration [cm/yr]');
    xlabel('Model time [Ma]'); title('Westward hinge motion'); set(gcf,'color','w');ylim([-2 4]);xlim([0 44]);
    
    total_retreat= indexmin_topo_new(closestIndexfull)-indexmin_topo_new;
    total_retreat(1:closestIndexfull-1)=0;
    
    cd(path_shortening)
    %Plot all the data
    load('oncken2012.mat');
    cummulative_shortening_oncken_2012=cumsum(oncken2012(:,2));
    
    %         shortening_oncken_2012_interp=interp1(oncken2012(:,1),oncken2012(:,2),0:0.05:50);
    %     shortening_oncken_2012_interp(isnan(shortening_oncken_2012_interp)) = [];cummulative_shortening_oncken_2012=cumsum(shortening_oncken_2012_interp);
    
    
    figure(15);clf;plot(flip(time./1e6),total_retreat);hold on;plot(oncken2012(:,1),cummulative_shortening_oncken_2012);set ( gca, 'xdir', 'reverse' );
%     title('Cumulative model trench retreat', 'Cumulative Oncken 2012 shortening')
    %     figure(17);clf;plot(oncken2012(:,1),flip(oncken2012(:,2)));
    
    figure(17);plot(oncken2012(:,1),oncken2012(:,2));
end

if any(Output_numbers==[1,2,8])
    figure(2); clf;
    if(isnumeric(plot_topography_time))
        index_map=plot_topography_time./1e6/(dt./1e6);
        indexmin_model=find(Elevationmap(index_map,:)==min(Elevationmap(index_map,:)));
        plot(topogreatcircle(:,3)-topogreatcircle(indexmin_topo,3),topocorrected/1e3,'LineWidth',line_width);hold on ;
        plot(Vec_lenght-Vec_lenght(indexmin_model),(Elevationmap(index_map,:)-bathymetry)/1e3,'LineWidth',line_width);
        legend('Real topography corrected from water','Model');
        xlabel('Distance from the trench [km]');
        ylabel('Altitude [km]'); title('Present day topography'); set(gcf,'color','w');%set(gca, 'color', 'none');
    elseif(ischar(plot_topography_time))
        % plot(topogreatcircle(:,3)-topogreatcircle(indexmin_topo,3),topogreatcircle(:,end),'LineWidth',line_width); hold on ;
        plot(topogreatcircle(:,3)-topogreatcircle(indexmin_topo,3),topocorrected,'LineWidth',line_width);hold on ;
        plot(Vec_lenght-Vec_lenght(indexmin_model(1)),Elevationmap(end,:)-bathymetry,'LineWidth',line_width);
        legend('Real topography corrected from water','Model');
        xlabel('Distance from the trench [km]');
        ylabel('Altitude [km]'); title('Present day topography'); set(gcf,'color','w');%set(gca, 'color', 'none');
    end
end



if any(Output_numbers==[10])
    cd(path_model);
    filenew = dir ('topolayer.*');
    if(isempty(filenew))
        return;
    else
        filenames_new = {filenew.name};
        filenames_new = natsortfiles(filenames_new);
        B=numel(filenames_new);
        
        cut_depth=700e3;
        
        incertitude_interval=5e3;
        % lets try to get an average dip with interval of depth depending on x interval
        % interval=50e3;
        % interval_depth_extract=interval:interval:cut_depth;
        interval_depth_extract=[35e3;75e3;150e3;300e3];
        % interval_depth_extract(end)=interval_depth_extract(end)-incertitude_interval;
        interval_depth_extract(end)=interval_depth_extract(end);
        min_interval_depth_extract=interval_depth_extract-incertitude_interval;
        max_interval_depth_extract=interval_depth_extract+incertitude_interval;
        % closestindex=1;
       
        
        clear z;
        z=0;
        for i=0:resampling_topo:B
            file_index=i;
            if (i==0)
                file_index=1;
            end
            z=z+1;
            pathelev = sprintf('%s',pwd,'/',filenames_new{file_index});
            fprintf('%i/%i\n', i,B)
            
            
            %% Setup the Import Options and import the data
            opts = delimitedTextImportOptions("NumVariables", 2);
            
            % Specify range and delimiter
            opts.DataLines = [2, Inf];
            opts.Delimiter = " ";
            
            % Specify column names and types
            opts.VariableNames = ["VarName1", "x"];
            opts.VariableTypes = ["double", "double"];
            
            % Specify file level properties
            opts.ExtraColumnsRule = "ignore";
            opts.EmptyLineRule = "read";
            opts.ConsecutiveDelimitersRule = "join";
            opts.LeadingDelimitersRule = "ignore";
            
            % Import the data
            topolayer = readtable(pathelev, opts);
            
            %% Convert to output type
            topolayer = table2array(topolayer);
            
            maxdepth=max(max(topolayer(:,2)));
            
            [val_sorted index_sorted]= sortrows(topolayer,[1 2]);
            [Xunique_sorted, Xunique_index] = unique(val_sorted(:,1),'last');
            elevations1=val_sorted(Xunique_index,2);
            interp_vec=0:1e3:modellenght*1e3;
            elevations2=interp1(Xunique_sorted,elevations1,interp_vec);
            Elevationmap2(z,:)=modelheight.*1e3-elevations2;
            
            elevation2_depth=modelheight.*1e3-elevations2;
            elevations2_cut_index = find(elevation2_depth<cut_depth);
            elevations2_cut = elevations2(elevations2_cut_index);
            index_nan=find(isnan(elevations2_cut));
            smooth_topo=smooth(elevations2_cut,smooth_dip)';
            smooth_topo(index_nan)=0;
            
            
            % elevations2_cut{i,:}=elevations2(elevations2>250e3)
%             dipsmooth=atand(abs((diff(modelheight.*1e3-smooth_topo))./1e3));
              % whithout abs allows to see that the slab has negative dipping during flat slab
              dipsmooth=atand((diff(modelheight.*1e3-smooth_topo)./1e3));
              
            % first and last values should be nan otherwise the smoothing send back a wrong dip
            % dipsmooth(1:5)=nan;dipsmooth(:,end)=nan;
            
            % we can compare to the dip without smoothing
            % dip=atand(abs((diff(modelheight.*1e3-elevations2))./(diff(interp_vec))));
            % dip(1:5)=nan;
            
%             dipmap should be the same length
            dipsmooth(numel(dipsmooth):modellenght)=nan;
            dipmap(z,:)=dipsmooth;
            
            
            % statistic for depth intervals
            for j=1:numel(interval_depth_extract)
                mean_diff_vec=abs(elevation2_depth-interval_depth_extract(j));
                [mean_closestvalue(j), mean_closestindex(j)]=find( mean_diff_vec == min(mean_diff_vec),1);
                min_diff_vec=abs(elevation2_depth-(interval_depth_extract(j)+incertitude_interval));
                [min_closestvalue(j), min_closestindex(j)]=find( min_diff_vec == min(min_diff_vec),1);
                max_diff_vec=abs(elevation2_depth-(interval_depth_extract(j)-incertitude_interval));
                [max_closestvalue(j), max_closestindex(j)]=find(max_diff_vec == min(max_diff_vec),1);
                
                %     if the interval we choose is greater than the data then resize it to fit the data number.
                if(min_closestindex(j)>numel(dipsmooth))
                    min_closestindex(j)=numel(dipsmooth);
                end
                
            end
            
            dipsmooth(isnan(dipsmooth))=0;
            for k=1:numel(mean_closestindex)
                nsample = min_closestindex(k)-max_closestindex(k);
                average_dipsmooth(k)=mean(dipsmooth(max_closestindex(k): min_closestindex(k)));
                ySEM = std(dipsmooth(max_closestindex(k): min_closestindex(k)))/sqrt(nsample);                              % Compute ‘Standard Error Of The Mean’ Of All Experiments At Each Value Of ‘x’
                CI95 = tinv([0.025 0.975], nsample-1);                    % Calculate 95% Probability Intervals Of t-Distribution
                yCI95{k} = bsxfun(@times, ySEM, CI95(:))';              % Calculate 95% Confidence Intervals Of All Experiments At Each Value Of ‘x’
            end
            average_dipsmooth_time(z,:)=average_dipsmooth;
            conf_interval(z,:)= yCI95;
        end
        %The timestep depends on what the user has decided in the Aspect prm file for the postprocess
        line_width =2;
        time=0:final_model_time/(B/resampling_topo):final_model_time
        bathymetry=3e3;
        
        figure(8);clf;
        surf(1:size(Elevationmap2,2),time,Elevationmap2./1e3);shading interp;c=colorbar; colormap(flipud(summer));ylabel('Time[Ma]'),xlabel('Length[km]');zlabel('Depth[km]');xlim([0 modellenght]);set(gcf,'color','w');set(gca, 'ZDir','reverse');%ylim([0 model_time_scaled])
        c.Label.String= "Slab surface depth [km]";view(2);title('Slab surface depth over time');%set(gca, 'color', 'none');grid off;set(gca,'XColor', 'none','YColor','none','ZColor','none');
        % % FaceLighting = 'gour
        
        figure(9);clf;
        surf(1:size(dipmap,2),time,dipmap);shading interp;c=colorbar; colormap(flipud(summer));ylabel('Time[Ma]'),xlabel('Length[km]');zlabel('Dip[°]');xlim([0 modellenght]);set(gcf,'color','w');set(gca, 'ZDir','reverse');%ylim([0 model_time_scaled])
        c.Label.String= "Dip [°]";view(2);title('Slab dip over time')%set(gca, 'color', 'none');grid off;set(gca,'XColor', 'none','YColor','none','ZColor','none');
        % % FaceLighting = 'gour
        
        % plot average dip at depth given depending on interval of uncertainties
        for k=1:size(average_dipsmooth_time,1)
            for i=1:size(average_dipsmooth_time,2)
%                 interval can be diproportionated if there are some nan
%                 value before calculation due to the absence of slab. So
%                 we chack if it is the case so they should not be calculated 
                if any(conf_interval{k,i}<-45)
                    conf_interval{k,i}(:)=nan;
                end
                
                full_interval_min(k,i)=conf_interval{k,i}(1)+average_dipsmooth_time(k,i);
                full_interval_max(k,i)=conf_interval{k,i}(2)+average_dipsmooth_time(k,i);
                
            end
        end
        
        figure(10);clf;
        color_depth= [204, 0, 153;
            0, 0, 204;
            102, 102, 153 ;
            255, 153, 0];
        X=[time,fliplr(time)];
        for i=1:size(average_dipsmooth_time,2)
        
%             we need to make a correction when the slab had not reached
%             the interval of depth yet
            index_nan=(~isnan(full_interval_min(:,i)));
            X_index_nan=[index_nan', fliplr(index_nan')];
            X_time=X(X_index_nan);
            %               plot(time,full_interval_min(:,i));
            %              plot(time,full_interval_max(:,i));
           
            
           full_interval_min_corrected=full_interval_min(index_nan,i) ;
           full_interval_max_corrected=full_interval_max(index_nan,i) ;
           
            full_interval_min_corrected(:)', flip(full_interval_max_corrected(:))';
            Y=[full_interval_min_corrected(:)', flip(full_interval_max_corrected(:))'];
            % fill(X,Y,,'Color',color_depth(i,:));
            patch(X_time,Y,1,'FaceColor',color_depth(i,:)./255,'EdgeColor','none');alpha(.2);hold on; % make patch transparent
            
%             Additionally we don'twant to plot the dip if the full
%             interval of depth is not reache by the slab
            index_no_slab=(index_nan==0); 
            average_dipsmooth_time(index_no_slab,i)=nan;
            
        end
        
        for k=1:size(average_dipsmooth_time,1)
            for i=1:size(average_dipsmooth_time,2)
                plot(time(k),average_dipsmooth_time(k,i),'ko-');
%                 plot(time(k),average_dipsmooth_time(k,i));
            end
        end
        legend('Dip 35±5km depth','Dip 75±5km depth','Dip 150±5km depth','Dip 300±5km depth');xlabel('Time[Ma]'),ylabel('Dip[°]');set(gcf,'color','w');title('Average dip in depth over time')%ylim([0 model_time_scaled])
        %set(gca, 'color', 'none');grid off;set(gca,'XColor', 'none','YColor','none','ZColor','none');
        % % FaceLighting = 'gour
         
        figure(20);clf;
       for i=1:size(average_dipsmooth_time,2)
       plot(time,average_dipsmooth_time(:,i));hold on
        end
        legend('Dip 35±5km depth','Dip 75±5km depth','Dip 150±5km depth','Dip 300±5km depth');xlabel('Time[Ma]'),ylabel('Dip[°]');set(gcf,'color','w');title('Average dip in depth over time');

        %comparison with real slab surface
        
        
        opts = delimitedTextImportOptions("NumVariables", 3);
        % Specify range and delimiter
        opts.DataLines = [1, Inf];
        opts.Delimiter = ",";
        % Specify column names and types
        opts.VariableNames = ["VarName1", "VarName2", "NaN"];
        opts.VariableTypes = ["double", "double", "double"];
        % Specify file level properties
        opts.ExtraColumnsRule = "ignore";
        opts.EmptyLineRule = "read";
        % Import the data for depth
        surfaceslab = readtable("/home/ponsm/Nextcloud/phd-central-andes-shortening/cookbookmika/dip_slab/sam_slab2_dep_02.23.18.xyz", opts);
        surfaceslab = table2array(surfaceslab);
        data_lat=surfaceslab(surfaceslab(:,2)==-21,:);
        data_lat=data_lat(~isnan(data_lat(:,3)),:);
        Xlon=deg2km(data_lat(:,1))-deg2km(data_lat(1,1));
        % Import the data for  uncertainties
        surfaceslab_uncertainties = readtable("/home/ponsm/Nextcloud/phd-central-andes-shortening/cookbookmika/dip_slab/sam_slab2_unc_02.23.18.xyz", opts);
        surfaceslab_uncertainties = table2array(surfaceslab_uncertainties);
        data_lat_uncertainties=surfaceslab_uncertainties(surfaceslab_uncertainties(:,2)==-21,:);
        data_lat_uncertainties=data_lat_uncertainties(~isnan(data_lat_uncertainties(:,3)),:);
        % Xlon_uncertainties=deg2km(data_lat_uncertainties(:,1))-deg2km(data_lat_uncertainties(1,1));
        % make curves for min and max
        max_depth_slab=abs(data_lat(:,3))+data_lat_uncertainties(:,3)./2;
        min_depth_slab=abs(data_lat(:,3))-data_lat_uncertainties(:,3)./2;
        X2=[Xlon',fliplr(Xlon')];
        Y2=[min_depth_slab(:)', flip(max_depth_slab(:)')];
        
        Vec_lenght= 1:size(Elevationmap2,2);
        
            %To georeference then lets use the minimum elevation, the trench assuming 8km depth.
        indexmin_model=find(abs(Elevationmap2(end,:) - 8e3)==min(abs(Elevationmap2(end,:) - 8e3)));
        
        figure(11);clf;
        patch(X2,Y2,1,'FaceColor',color_depth(1,:)./255,'EdgeColor','none');alpha(.2); hold on;% make patch transparent
        plot(Xlon,abs(data_lat(:,3)));
        % plot(Xlon_uncertainties,max_depth_slab,Xlon_uncertainties,min_depth_slab)
        plot (Vec_lenght-Vec_lenght(indexmin_model(1)),Elevationmap2(end,:)./1e3);
        legend('Slab2 uncertainties','Slab2 mean','Slab (Model)');xlabel('Distance from trench[km]');ylabel('Depth[km]');xlim([0 max(Xlon)]);set(gcf,'color','w');set(gca, 'YDir','reverse');hold on;%ylim([0 model_time_scaled])
        title('Present day slab surface')
        
        
        Vec_dip_smooth=1:numel(dipsmooth);
        % Import the data dip data from slab2
        dip_slab = readtable("/home/ponsm/Nextcloud/phd-central-andes-shortening/cookbookmika/dip_slab/sam_slab2_dip_02.23.18.xyz", opts);
        dip_slab  = table2array(dip_slab );
        data_lat_dip_slab =dip_slab (dip_slab (:,2)==-21,:);
        data_lat_dip_slab =data_lat_dip_slab (~isnan(data_lat_dip_slab (:,3)),:);
        figure(12);clf;
        plot(Xlon,data_lat_dip_slab(:,3)); hold on;
        plot(Vec_dip_smooth-Vec_lenght(indexmin_model(1)),dipsmooth);set(gca, 'YDir','reverse');xlabel('Distance from trench[km]');ylabel('Dip[°]');xlim([0 max(Xlon)]);set(gcf,'color','w');title('Present day slab dip');
        legend('Slab2','Slab (Model)');
        
    end
end



ind_ref_start = find(min(abs(Data(:,2) - 6.5e6)) == abs(Data(:,2) - 6.5e6),1,'last');
ind_ref_stop = find(min(abs(Data(:,2) - 38e6)) == abs(Data(:,2) - 38e6),1,'last');
Time_account = Data(ind_ref_start:ind_ref_stop,2); 
vel_account = Data(ind_ref_start:ind_ref_stop,data_number);

ind_ref_start_Vt = find(min(abs(time - 6.5e6)) == abs(time - 6.5e6),1,'last');
ind_ref_stop_Vt = find(min(abs(time - 38e6)) == abs(time - 38e6),1,'last');
Vt_time_account = time(ind_ref_start_Vt:ind_ref_stop_Vt); 
Vt_account = trench_rate_smooth(ind_ref_start_Vt:ind_ref_stop_Vt);

% Vd=Data(ind_ref,shortening_stat)-Data(:,shortening_stat);
% Vdf=Data(:,underthrusting_stat)-Data(ind_ref,underthrusting_stat); 
% Vdt = Vd+Vdf;
% Vd_account= Vd(ind_ref_start:ind_ref_stop);
% Vdf_account = Vdf(ind_ref_start:ind_ref_stop);
% Vdt_account = Vdt(ind_ref_start:ind_ref_stop);
% dtime_account=diff(Time_account); 
% Vdt_rate=diff(Vdt_account.*1e3)./dtime_account
% Vdt_rate_smooth=smooth(Vdt_rate,smoothing_stats_interval);

ind_ref_start_Vdt = find(min(abs(time_vec - 6.5e6)) == abs(time_vec - 6.5e6),1,'last');
ind_ref_stop_Vdt = find(min(abs(time_vec - 38e6)) == abs(time_vec - 38e6),1,'last');
Vdt_time_account = time_vec(ind_ref_start_Vdt:ind_ref_stop_Vdt); 
total_shortening_rate_smooth_adj = [0 total_shortening_rate_smooth'];
Vdt_account = total_shortening_rate_smooth_adj(ind_ref_start_Vdt:ind_ref_stop_Vdt);


tt=6.5e6:0.1e6:38e6; 
Vs_interp = interp1(Time_account,vel_account,tt);
Vt_interp = interp1(Vt_time_account,Vt_account,tt);
% Vdt_rate_smooth= interp1(Time_account(2:end),Vdt_rate_smooth,tt);
Vdt_interp=interp1(Vdt_time_account,Vdt_account,tt);

clear color_phases str_symbol Vs_grp Vt_grp Vdt_grp symb color_top phases  color_phases  ind_colors
% phases =[6.5e6;10e6;12e6;15e6;20e6;26e6;29e6;31e6;35e6;38e6];
 phases =[6.5e6;11e6;20e6;29e6;38e6];

for i=1:numel(phases)
    ind_colors(i) = find(min(abs(tt - phases(i))) == abs(tt - phases(i)),1,'last');
    if (ind_colors(i)==0)
        ind_colors(i)=1;
    end
end
% symb=[{'o'} {'p'} {'^'} {'s'} {'h'}];
for i=2:numel(ind_colors)
    color_phases(ind_colors(i-1):ind_colors(i))=i-1;
    Vs_grp{i-1}=Vs_interp(ind_colors(i-1):ind_colors(i)).*1e2;
    Vt_grp{i-1}=Vt_interp(ind_colors(i-1):ind_colors(i)); 
    Vdt_grp{i-1}=Vdt_interp(ind_colors(i-1):ind_colors(i)); 
end

% symb=['h'; 'p'; '^' ;'s';'o';'.';'*';'x';'+'];
% color_top =['k';'k';'k';'k';'k';'k';'k';'k';'k'];
symb=['o'; 'p'; '^' ;'s';'h'];
color_top =['k';'k';'k';'k';'k'];
figure (16); clf; 
for i=1:max(color_phases)
    hold on;
% scatter(Vs_grp{i},Vt_grp{i},50,Vdt_grp{i},'filled',symb(i));xlabel('VSub[cm/yr]');ylabel('Vt[cm/yr]');set(gcf,'color','w');title('Vt/Vsub');colorbar;
% ultime plot colored filled edeg 
scatter(Vs_grp{i},Vt_grp{i},60,Vdt_grp{i},'filled',symb(i),'MarkerEdgeColor',color_top(i),'LineWidth',1);xlabel('VSub[cm/yr]');xlabel('VSub[cm/yr]');ylabel('Vt[cm/yr]');set(gcf,'color','w');title('Vt/Vsub');c=colorbar;
end
crameri('lajolla');
c.Label.String= "Shortening rate [mm/yr]"
legend('Phase 1','Phase 2','Phase 3','Phase 4');


figure(17); 
plot(Vdt_interp)

end
