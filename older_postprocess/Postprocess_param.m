% Postprocesses 2 

%% Collect_data
clc;clear all;close all; 
% As script that collect all data needed from the subduction models
addpath(pwd)

%% Output options
% 1-Output topography
% 2-Output topography at X Ma
% 3-Output topography profil evolution
% 4-Output Migration trench rate
% 5-Output Oceanic plate velocity
% 6-Output shortening rate and shortening
% 7-Output Value for cumulated shortening before X Ma
% 8-plot All
% 9-plot slab dip

set(0,'defaultfigureposition',[500 500 400 300]')
use_old_parameters=false; 

%% Define the path

% cd '/home/ponsm/Desktop/model_g/model/C_model_28/571_549_M5b_vel_OP4cm/'
cd '/Users/ponsm/Desktop/modelblogin/model/globalscale/anulus2d/mobility_function/output'

path_model=pwd;

    %% parameters
% % Statistic number for the velocity or other statistics to plot
% Statistic number for the shortening
% 
data_number=81;  
shortening_stat= 88; 
underthrusting_stat=89;


color_time =[6.5e6;11e6;20e6;29e6];
colors= [250, 235, 245;
    242, 242, 237;
    222, 235, 250 ;
    228, 237, 225;
    249, 250, 235];


shade_factor =0;
colors= (colors* (1 - shade_factor))/255;

%Give a list of output number
Output_numbers = [8]


%Default values

%Model scale time
% model_time_scaled=44e6;
model_time_scaled=38e6;
% 1-dt topography
dt = 20000; 
add_OP_vel=0; %m/yr it should be 0 if we only want the abs motion

% 2- Time in Ma or a 'end' statement for the last
plot_topography_time=37e6;
% plot_topography_time='end';%default should be


% 6 - Initial time of reference for the shortening calculation
initial_state_time = 6.5e6
% 7-Total cumulative shortening 
cumultime=44e6;

%Plugin topograpgraphy can be smoothed too, resolution = outputed topography time interval from Prm *resampling_topo * smoothing interval
% if resampling topo is 1 all the files will be taken but the computation will be longer and ll require more memory
final_model_time= 40e6;
resampling_topo=1;
smoothing_interval = 25; 
smooth_dip = 50; 
%Resample stats to speedup the smoothing we can resample every 500e3 yr time_resolution = resampling * smoothing interval
resampling_resolution=500e3; smoothing_stats_interval = 10; 

path_shortening='/home/ponsm/Nextcloud/phd-central-andes-shortening/cookbookmika/shortenning'
%% Scheme

get_allstats(path_model,Output_numbers,cumultime,model_time_scaled,initial_state_time,path_shortening,data_number,shortening_stat,underthrusting_stat,color_time,colors,resampling_resolution,smoothing_stats_interval,plot_topography_time,dt,smoothing_interval,smooth_dip,resampling_topo,final_model_time,add_OP_vel)