clc ; close all; clear all;
addpath(genpath('./Individual_scripts/'))

% path_models = {'/Users/ponsm/Desktop/modelblogin/model/globalscale/Potsprocess/sphere3d/P01a_Pangea_1GPa_Mantle_C40MPa_LR/'}
path_models = {'/Users/ponsm/Desktop/Model_sphere/Models_HLRN/P01a_Pangea_1GPa_Mantle_C40MPa_LR/'}
%     '/Users/ponsm/Desktop/Model_sphere/Models_HLRN/R01f_Rodinia_2GPa_Mantle_C10MPa_f005_LR/'}
%     '/Users/ponsm/Desktop/Model_sphere/Models_HLRN/P01a_Pangea_1GPa_Mantle_C40MPa_LR/'}

% Extracting model titles
model_titles = cell(size(path_models));
for t = 1:numel(path_models)
% Split the path into parts
parts = strsplit(path_models{t}, '/');
% Extract the last part
model_titles{t} = parts{end-1};
end

%Load colors
% colors_topography = load('colors.mat');

for m = 1:numel(path_models)
    path_model = cell2mat(path_models(m));
    make_video='true';
    resample_dynamic_topography = 10;

    %     statistic_parameters = {'Time'};
    %     % Get the data and corresponding statistic numbers and indices
    %     [data_stats, stats_number, header, stat_indices] = get_statistics(path_model, statistic_parameters);
    %
    %     % Statistics
    %     time_vec = data_stats(:, 2)./1e6;
    %     time=num2str(time_vec, '%05.2f Myr');
    %
    %     % find how many zeros for the initial refinement steps
    %     numbers_of_IRS=nnz(~time_vec);

    file_surface = dir(fullfile(path_model, 'surface_*'));
    filenames_surface = {file_surface.name};
    num_files_surface = numel(filenames_surface);

%     file_lithosphere = dir(fullfile(path_model, 'lithosphere_*'));
%     filenames_lithosphere = {file_lithosphere.name};
%     num_files_lithosphere = numel(filenames_lithosphere);
% 
%     file_depths = dir(fullfile(path_model, 'depths_*'));
%     filenames_depths = {file_depths.name};
%     num_files_depths = numel(filenames_depths);

    % Specify the desired frame size for the video
    %     frameSize = [2000, 1000];

    % Create the figure with the specified size
    %     figure('WindowState', 'maximized','Position', [100 100 frameSize]);
    %     scatterSize =70;

    if strcmp(make_video, 'true')
        % Specify the video filename
        videoFilename = [model_titles{m} '.mp4'];

        % Create a VideoWriter object
        videoWriter = VideoWriter(videoFilename, 'MPEG-4');
        videoWriter.FrameRate = 4; % Specify the frame rate (adjust as needed)
        open(videoWriter);
    end

%     count_ite = 0;
    Reference_mean_correction_steps=2;

    %     We want the initial position of the continent to be plotted on top of the maps
    surface_init = import_surface_sphere(fullfile(path_model, filenames_surface{Reference_mean_correction_steps}));
    % Sort the data in the desired order
    x_surface_init = surface_init.Points0;
    y_surface_init = surface_init.Points1;
    z_surface_init = surface_init.Points2;

    % Calculate the radial distance from the center of the Earth
    r_surface_init = sqrt(x_surface_init.^2 + y_surface_init.^2 + z_surface_init.^2);

    surface_topography_init = r_surface_init-6371e3;

    % Calculate the azimuthal angle
    theta_surface_init = atan2(y_surface_init, x_surface_init);
    [theta_sorted_init, theta_index_init] = sort(theta_surface_init);
    x_sorted_surface_init = x_surface_init(theta_index_init);
    y_sorted_surface_init = y_surface_init(theta_index_init);
    z_sorted_surface_init = z_surface_init(theta_index_init);
     surface_init_topography_sorted = surface_topography_init(theta_index_init);
    surface_init_continent_sorted = surface_init.continent(theta_index_init);
    % Convert Cartesian coordinates to longitude and latitude
    longitude_surface_init = atan2(y_sorted_surface_init, x_sorted_surface_init) * 180 / pi;
    latitude_surface_init = asin(z_sorted_surface_init ./ sqrt(x_sorted_surface_init.^2 + y_sorted_surface_init.^2 + z_sorted_surface_init.^2)) * 180 / pi;
    % Create a grid of longitude and latitude values
    [Xeq, Yeq] = meshgrid(-180:1:180, -90:1:90);
    [Xeq_lr, Yeq_lr] = meshgrid(-180:5:180, -90:5:90);
    continents_init = griddata(longitude_surface_init, latitude_surface_init, surface_init_continent_sorted, Xeq, Yeq);

    mean_surface_init_topography = mean(surface_init_topography_sorted);
%     for i =99
            for i = 1:num_files_surface

        count_ite = i;
        fprintf('Processing surface file %d of %d\n', i, num_files_surface)
        file_index = i;

        % Import the data
        surface = import_surface_sphere(fullfile(path_model, filenames_surface{file_index}));
%         lithosphere = import_lithosphere_sphere(fullfile(path_model, filenames_lithosphere{file_index}));
%         depths = import_depths_sphere(fullfile(path_model, filenames_depths{file_index}));

        %         filename_extract = filenames{file_index};
        %         parts = split(filename_extract, '_');
        %         number_step = str2double(parts{2});
        time_str = sprintf('Time: %.2f Myr', surface.Time(1)./1e6);
        % %05.2f Myr

        % Define the downsampling factor
        downsample_factor = 1;

        % Sort the data in the desired order
        x_surface = surface.Points0;
        y_surface = surface.Points1;
        z_surface = surface.Points2;

%         x_lithosphere = lithosphere.Points0;
%         y_lithosphere = lithosphere.Points1;
%         z_lithosphere = lithosphere.Points2;
% 
%         x_depths = depths.Points0;
%         y_depths = depths.Points1;
%         z_depths = depths.Points2;

        % Calculate the radial distance from the center of the Earth
        r_surface = sqrt(x_surface.^2 + y_surface.^2 + z_surface.^2);

        surface_topography = r_surface-6371e3;

        % Calculate the azimuthal angle
        theta_surface = atan2(y_surface, x_surface);

        % Sort the data based on the azimuthal angle
        [theta_sorted, theta_index] = sort(theta_surface);
        x_sorted_surface = x_surface(theta_index);
        y_sorted_surface = y_surface(theta_index);
        z_sorted_surface = z_surface(theta_index);
        r_surface_sorted = r_surface(theta_index);
        surface_topography_sorted = surface_topography(theta_index);
        surface_strain_rate_sorted = surface.strain_rate(theta_index);
        surface_continent_sorted = surface.continent(theta_index);
        surface_v_phi_sorted = surface.v_phi(theta_index);
        surface_v_theta_sorted = surface.v_theta(theta_index);
        surface_v_r_sorted = surface.v_r(theta_index);
        surface_velocity0_sorted = surface.velocity0(theta_index);
        surface_velocity1_sorted = surface.velocity1(theta_index);
        surface_velocity2_sorted = surface.velocity2(theta_index);

%         theta_lithosphere = atan2(y_lithosphere, x_lithosphere);
%         [theta_sorted_lithosphere, theta_index_lithosphere] = sort(theta_lithosphere);
%         x_sorted_lithosphere = x_lithosphere(theta_index_lithosphere);
%         y_sorted_lithosphere = y_lithosphere(theta_index_lithosphere);
%         z_sorted_lithosphere = z_lithosphere(theta_index_lithosphere);
%         divergence_lithosphere_sorted = lithosphere.Divergence(theta_index_lithosphere);
%         oceanic_age_lithosphere_sorted = lithosphere.oceanic_age(theta_index_lithosphere);
% 
%         theta_depths = atan2(y_depths, x_depths);
%         [theta_sorted_depths, theta_index_depths] = sort(theta_depths);
%         x_sorted_depths = x_depths(theta_index_depths);
%         y_sorted_depths = y_depths(theta_index_depths);
%         z_sorted_depths = z_depths(theta_index_depths);
%         depths_sorted_depths = depths.depth(theta_index_depths);
%         unique_depths = unique(depths_sorted_depths);
%         depths_nonadiabatic_temperature_sorted = depths.nonadiabatic_temperature(theta_index_depths);
%         depths_v_phi_sorted = depths.v_phi(theta_index_depths);
%         depths_v_theta_sorted = depths.v_theta(theta_index_depths);
%         depths_v_r_sorted = depths.v_r(theta_index_depths);
%         depths_velocity0_sorted = depths.velocity0(theta_index_depths);
%         depths_velocity1_sorted = depths.velocity1(theta_index_depths);
%         depths_velocity2_sorted = depths.velocity2(theta_index_depths);

%         for u = 1:numel(unique_depths)
%             idx_depths = find(depths_sorted_depths == unique_depths(u));
%             x_sorted_depths_map=x_sorted_depths(idx_depths) ;
%             y_sorted_depths_map=y_sorted_depths(idx_depths);
%             z_sorted_depths_map=z_sorted_depths(idx_depths);
%             x_depths_downsampled{u} = x_sorted_depths_map(1:downsample_factor:end);
%             y_depths_downsampled{u} = y_sorted_depths_map(1:downsample_factor:end);
%             z_depths_downsampled{u} = z_sorted_depths_map(1:downsample_factor:end);
%             non_adiabT_map{u}=depths_nonadiabatic_temperature_sorted(idx_depths);
%             v_theta_map{u}=depths_v_theta_sorted(idx_depths);
%             v_phi_map{u}=depths_v_phi_sorted(idx_depths);
%             v_r_map{u}=depths_v_r_sorted(idx_depths);
%             velocity0_map{u}=depths_velocity0_sorted(idx_depths);
%             velocity1_map{u}=depths_velocity1_sorted(idx_depths);
%             velocity2_map{u}=depths_velocity2_sorted(idx_depths);
%         end

        %         final_surface_topography_sorted = surface_topography_sorted - ((surface_topography_sorted < 2400) .* (1000 / 3300) .* abs(surface_topography_sorted - 2400)) - 2400;
        final_surface_topography_sorted=surface_topography_sorted;

% 
%         x_surface_downsampled = x_sorted_surface(1:downsample_factor:end);
%         y_surface_downsampled = y_sorted_surface(1:downsample_factor:end);
%         z_surface_downsampled = z_sorted_surface(1:downsample_factor:end);
%         surface_topography_surface_downsampled = final_surface_topography_sorted(1:downsample_factor:end);

%         x_lithosphere_downsampled = x_sorted_lithosphere(1:downsample_factor:end);
%         y_lithosphere_downsampled = y_sorted_lithosphere(1:downsample_factor:end);
%         z_lithosphere_downsampled = z_sorted_lithosphere(1:downsample_factor:end);


        min_surface_topography = min(final_surface_topography_sorted);
        max_surface_topography = max(final_surface_topography_sorted);
        mean_surface_topography(count_ite) = mean(final_surface_topography_sorted);

        %        Correction for mean topography subsiding
        if(count_ite<=Reference_mean_correction_steps)
            diff_mean_surface_topography(count_ite)=0;
            disp('test1')
        elseif (mean_surface_topography(count_ite)<mean_surface_topography(count_ite-1))
            diff_mean_surface_topography(count_ite)=mean_surface_init_topography-mean_surface_topography(count_ite);
            disp('test2')
        else
            diff_mean_surface_topography(count_ite)=0;
            disp('test3')
        end


%         % Convert Cartesian coordinates to longitude and latitude
        longitude_surface_downsampled = atan2(y_sorted_surface, x_sorted_surface) * 180 / pi;
        latitude_surface_downsampled = asin(z_sorted_surface ./ sqrt(x_sorted_surface.^2 + y_sorted_surface.^2 + z_sorted_surface.^2)) * 180 / pi;

%         longitude_lithosphere_downsampled = atan2(y_lithosphere_downsampled, x_lithosphere_downsampled) * 180 / pi;
%         latitude_lithosphere_downsampled = asin(z_lithosphere_downsampled ./ sqrt(x_lithosphere_downsampled.^2 + y_lithosphere_downsampled.^2 + z_lithosphere_downsampled.^2)) * 180 / pi;
% 


        % Plot the data on a surface
        %     h2=subplot(1, 2, 2);
        % Create a grid of longitude and latitude values
%         [Xeq, Yeq] = meshgrid(-180:1:180, -90:1:90);
        [Xeq_lr, Yeq_lr] = meshgrid(-180:5:180, -90:5:90);
        % Reshape the surface topography data to match the grid size
        Topo = griddata(longitude_surface_downsampled, latitude_surface_downsampled, final_surface_topography_sorted, Xeq, Yeq);
        Topo = Topo+diff_mean_surface_topography(count_ite);
%         V_phi_vec_surface = griddata(longitude_surface_downsampled, latitude_surface_downsampled,surface_v_phi_sorted, Xeq_lr, Yeq_lr);
%         V_theta_vec_surface = griddata(longitude_surface_downsampled, latitude_surface_downsampled,surface_v_theta_sorted, Xeq_lr, Yeq_lr);
%         V_phi_surface = griddata(longitude_surface_downsampled, latitude_surface_downsampled,surface_v_phi_sorted, Xeq, Yeq);
%         V_theta_surface = griddata(longitude_surface_downsampled, latitude_surface_downsampled,surface_v_theta_sorted, Xeq, Yeq);
%         V_r_surface = griddata(longitude_surface_downsampled, latitude_surface_downsampled,surface_v_r_sorted, Xeq, Yeq);
%         V_magnitude_surface = sqrt(V_r_surface.^2+V_phi_surface.^2+V_theta_surface.^2);
% 

%         strain_rate = griddata(longitude_surface_downsampled, latitude_surface_downsampled, surface_strain_rate_sorted, Xeq, Yeq);
        continents = griddata(longitude_surface_downsampled, latitude_surface_downsampled, surface_continent_sorted, Xeq, Yeq);

%         for uu = 1:numel(unique_depths)
%             longitude_depths_downsampled = atan2(y_depths_downsampled{uu}, x_depths_downsampled{uu}) * 180 / pi;
%             latitude_depths_downsampled = asin(z_depths_downsampled{uu} ./ sqrt(x_depths_downsampled{uu}.^2 + y_depths_downsampled{uu}.^2 + z_depths_downsampled{uu}.^2)) * 180 / pi;
%             non_adiabT_depths{uu}=griddata(longitude_depths_downsampled, latitude_depths_downsampled, non_adiabT_map{uu}, Xeq, Yeq);
%             V_theta_vec_depths{uu}=griddata(longitude_depths_downsampled, latitude_depths_downsampled, v_theta_map{uu}, Xeq_lr, Yeq_lr);
%             V_phi_vec_depths{uu}=griddata(longitude_depths_downsampled, latitude_depths_downsampled, v_phi_map{uu}, Xeq_lr, Yeq_lr);
%             V_theta_depths{uu}=griddata(longitude_depths_downsampled, latitude_depths_downsampled, v_theta_map{uu}, Xeq, Yeq);
%             V_phi_depths{uu}=griddata(longitude_depths_downsampled, latitude_depths_downsampled, v_phi_map{uu}, Xeq, Yeq);
%             V_r_depths{uu}=griddata(longitude_depths_downsampled, latitude_depths_downsampled, v_r_map{uu}, Xeq, Yeq);
%             Velocity0_depths{uu}=griddata(longitude_depths_downsampled, latitude_depths_downsampled, velocity0_map{uu}, Xeq, Yeq);
%             Velocity1_depths{uu}=griddata(longitude_depths_downsampled, latitude_depths_downsampled, velocity1_map{uu}, Xeq, Yeq);
%             Velocity2_depths{uu}=griddata(longitude_depths_downsampled, latitude_depths_downsampled, velocity2_map{uu}, Xeq, Yeq);
%             V_magnitude_depths{uu}=sqrt(V_r_depths{uu}.^2+V_phi_depths{uu}.^2+V_theta_depths{uu}.^2);
%         end

        %         We do it also for the lithosphere as the coordinate may slightly
        %         change when extracting from paraview which may lead to some
        %         visualisation artefact
        %         continents_lithosphere = griddata(longitude_lithosphere_downsampled, latitude_lithosphere_downsampled, surface_continent_sorted, Xeq, Yeq);

        % smoothing_data_coefficient = 0.01;
        % smooth_strain_rate = imgaussfilt(strain_rate,smoothing_data_coefficient);
        scaling_velocity_factor_surface=0.2;
        scaling_velocity_factor_depths=0.01;
        CMap_plume = [1 0 0 ;  1 1 1];
        CMap_subduction = [0 0 1 ;  1 1 1];
        CMap_boundaries = [ 0.0471 0.4392 0.2275;  1 1 1];
        CMap_continents = [ 0.4745 0.3647 0.3020;  1 1 1];
        CMap_continents_dark = [ 0 0 0;  1 1 1];

        % Track continents with a thrshold
        continents_threshold = continents>0.3;
        idx_is_nan_to_zero=isnan(continents_threshold);
        continents_threshold(idx_is_nan_to_zero)=0;
        oceans = 1 - continents_threshold;
        %         continents(continents == 0) = NaN;
        continents_reversed_position_for_plot = ind2rgb(oceans + 1, CMap_continents);
        continents_reversed_position_for_plot_dark = ind2rgb(oceans + 1, CMap_continents_dark);

        
%% Plot corrected topography
 %% Plot Strain rate with topography
%  Topo_corrected=Topo;
 smooth_Topo = imgaussfilt(Topo,0.3);

 Topo_corrected=smooth_Topo ; 
 idx_oceans=find(oceans == 1);
  idx_cont=find(oceans == 0);
Topo_corrected(idx_oceans) = Topo_corrected(idx_oceans)-1900 - (Topo_corrected(idx_oceans)<0).* (1000 / 3300).* abs(Topo_corrected(idx_oceans));

%%
% Topography_colors = [...
%       -11000          10           0         121;
%       -10500          26           0         137;
%       -10000          38           0         152;
%        -9500          27           3         166;
%        -9000          16           6         180;
%        -8500           5           9         193;
%        -8000           0          14         203;
%        -7500           0          22         210;
%        -7000           0          30         216;
%        -6500           0          39         223;
%        -6000          12          68         231;
%        -5500          26         102         240;
%        -5000          19         117         244;
%        -4500          14         133         249;
%        -4000          21         158         252;
%        -3500          30         178         255;
%        -3000          43         186         255;
%        -2500          55         193         255;
%        -2000          65         200         255;
%        -1500          79         210         255;
%        -1000          94         223         255;
%         -500         138         227         255;
%           0          16         123          48;
%          500         232         214         125;
%         1200         163          68           0;
%         1800         130          30          30;
%         2800         161         161         161;
%         3000         206         206         206;
%         3500         200         200         200;
%         4000         180         180         180;
%         4500         160         160         160;
%         5000          90          90          90];

Topography_colors_test = [...
      -11000          10           0         121;
      -10500          26           0         137;
      -10000          38           0         152;
       -9500          27           3         166;
       -9000          16           6         180;
       -8500           5           9         193;
       -8000           0          14         203;
       -7500           0          22         210;
       -7000           0          30         216;
       -6500           0          39         223;
       -6000          12          68         231;
       -5500          26         102         240;
       -5000          19         117         244;
       -4500          14         133         249;
       -4000          21         158         252;
       -3500          30         178         255;
       -3000          43         186         255;
       -2500          55         193         255;
       -2000          65         200         255;
       -1500          79         210         255;
       -1000          94         223         255;
        -500         138         227         255;
          0          16         123          48;
         500         232         214         125;
        1000	243	202	137;
        1500	230	184	88;
        2000	 130          30          30;
        2500	 164	144	25;
        3000	162	134	19;
        3500	217	166	39;
        4000         180         180         180;
        4500         160         160         160;
        5000          90          90          90];
Topography_colors= [
-10000	25.86909048	38.40738046	89.66856008;
-9922.48062	27.44599819	39.88142148	91.1577821;
-9844.96124	29.00165234	41.36714014	92.65152773;
-9767.44186	30.53267709	42.85279301	94.14790216;
-9689.922481	32.06443574	44.34525283	95.65420206;
-9612.403101	33.5933824	45.85193093	97.16322656;
-9534.883721	35.11492765	47.36630918	98.67700262;
-9457.364341	36.63467512	48.87898834	100.2014421;
-9379.844961	38.15305351	50.41073366	101.7262228;
-9302.325581	39.6755868	51.93703836	103.2610465;
-9224.806202	41.20425162	53.4768076	104.7991548;
-9147.286822	42.73667156	55.01997196	106.3436043;
-9069.767442	44.26445822	56.57986933	107.8933682;
-8992.248062	45.80759755	58.13515464	109.454383;
-8914.728682	47.3547645	59.6944798	111.0140692;
-8837.209302	48.9009794	61.27031397	112.5860258;
-8759.689922	50.4639815	62.84353741	114.1595321;
-8682.170543	52.02557375	64.43359984	115.7449032;
-8604.651163	53.59593505	66.01716865	117.3301098;
-8527.131783	55.17477584	67.61779347	118.926239;
-8449.612403	56.75770084	69.22030414	120.5243593;
-8372.093023	58.352945	70.82593868	122.1286205;
-8294.573643	59.94832221	72.43925951	123.7386641;
-8217.054264	61.54722666	74.05736124	125.3549849;
-8139.534884	63.16085063	75.67814853	126.9772727;
-8062.015504	64.77989664	77.30649046	128.6066049;
-7984.496124	66.40504958	78.94513156	130.2418728;
-7906.976744	68.03215891	80.58347018	131.8787453;
-7829.457364	69.66149551	82.23185475	133.5223269;
-7751.937984	71.30797848	83.88219998	135.1728378;
-7674.418605	72.9507805	85.54224329	136.8281533;
-7596.899225	74.60842708	87.20291864	138.4879724;
-7519.379845	76.26250832	88.87171689	140.1535986;
-7441.860465	77.92749974	90.54014343	141.8246726;
-7364.341085	79.60022124	92.2173887	143.5005091;
-7286.821705	81.27268139	93.90214878	145.1824369;
-7209.302326	82.95162261	95.58815158	146.8671102;
-7131.782946	84.63732134	97.28078945	148.5586829;
-7054.263566	86.33062713	98.97598387	150.2531929;
-6976.744186	88.02348837	100.6792415	151.9547875;
-6899.224806	89.72872246	102.3857426	153.6598805;
-6821.705426	91.43514366	104.0990277	155.3715395;
-6744.186047	93.14235708	105.8161415	157.0841142;
-6666.666667	94.8578028	107.5363731	158.8068724;
-6589.147287	96.58103731	109.2626711	160.5296206;
-6511.627907	98.30618676	110.9908351	162.2601515;
-6434.108527	100.0376081	112.7271124	163.9930216;
-6356.589147	101.7722738	114.4664961	165.7314045;
-6279.069767	103.5113187	116.2113002	167.4739473;
-6201.550388	105.257202	117.958515	169.2209352;
-6124.031008	107.0059086	119.7131642	170.9745817;
-6046.511628	108.7616857	121.4713965	172.7295415;
-5968.992248	110.5217734	123.23142	174.4891345;
-5891.472868	112.2813579	124.9995071	176.256304;
-5813.953488	114.0484743	126.7692616	178.0251393;
-5736.434109	115.8240203	128.5439332	179.7981745;
-5658.914729	117.599032	130.3233889	181.5776474;
-5581.395349	119.3783177	132.103696	183.3599455;
-5503.875969	121.1615558	133.8940269	185.1444473;
-5426.356589	122.9531573	135.6850908	186.9364601;
-5348.837209	124.7475961	137.4822581	188.7304354;
-5271.317829	126.5451047	139.2802411	190.5289398;
-5193.79845	128.3471753	141.0843178	192.3329913;
-5116.27907	130.1539474	142.8905096	194.1385354;
-5038.75969	131.9606251	144.7021256	195.9487896;
-4961.24031	133.7757564	146.520664	197.7607619;
-4883.72093	135.5952898	148.3409878	199.577233;
-4806.20155	137.4184821	150.1632952	201.3977949;
-4728.682171	139.2428092	151.9911744	203.2190157;
-4651.162791	141.073041	153.8218241	205.0403749;
-4573.643411	142.905038	155.6585962	206.8639082;
-4496.124031	144.7424507	157.4946845	208.6863489;
-4418.604651	146.5863341	159.336618	210.508047;
-4341.085271	148.4307655	161.1814332	212.3252468;
-4263.565891	150.2758764	163.030486	214.1356886;
-4186.046512	152.1260887	164.8778649	215.9377352;
-4108.527132	153.9771839	166.7302882	217.7295157;
-4031.007752	155.8320225	168.5807627	219.505955;
-3953.488372	157.6824207	170.4323296	221.2645972;
-3875.968992	159.5347209	172.2823851	222.9978174;
-3798.449612	161.383155	174.1259831	224.7034368;
-3720.930233	163.2257591	175.9679007	226.3752362;
-3643.410853	165.0611936	177.7997134	228.0061805;
-3565.891473	166.8901042	179.6231743	229.5931405;
-3488.372093	168.7033468	181.4360718	231.1274246;
-3410.852713	170.503925	183.2314043	232.6041235;
-3333.333333	172.2869798	185.0059703	234.0192385;
-3255.813953	174.0457578	186.7629205	235.3652793;
-3178.294574	175.7838688	188.4933243	236.6412801;
-3100.775194	177.4949963	190.2006295	237.8428629;
-3023.255814	179.1762854	191.8751605	238.9684434;
-2945.736434	180.8295974	193.5232864	240.0188721;
-2868.217054	182.448152	195.1382174	240.9915936;
-2790.697674	184.0379041	196.7222729	241.8920397;
-2713.178295	185.5963093	198.2756063	242.7195908;
-2635.658915	187.1212148	199.7990118	243.4812233;
-2558.139535	188.6186489	201.2943191	244.1800152;
-2480.620155	190.089072	202.7601717	244.8215382;
-2403.100775	191.5331496	204.2020394	245.4118788;
-2325.581395	192.9537832	205.6220581	245.9564448;
-2248.062016	194.3544853	207.023298	246.4607083;
-2170.542636	195.7392286	208.4055328	246.9304926;
-2093.023256	197.107653	209.7742238	247.3722437;
-2015.503876	198.4644934	211.1317963	247.7905871;
-1937.984496	199.8118678	212.4804592	248.1870906;
-1860.465116	201.1520987	213.8202871	248.5685335;
-1782.945736	202.4848093	215.1554684	248.9381546;
-1705.426357	203.8138183	216.4849527	249.2960428;
-1627.906977	205.1405744	217.8124921	249.6468403;
-1550.387597	206.4664015	219.1399129	249.9924875;
-1472.868217	207.7898686	220.4649302	250.3339618;
-1395.348837	209.1136301	221.7910519	250.6725228;
-1317.829457	210.4384207	223.1178042	251.0087832;
-1240.310078	211.7629145	224.4435858	251.3439545;
-1162.790698	213.0889955	225.773221	251.6770647;
-1085.271318	214.4170883	227.1030189	252.0113134;
-1007.751938	215.7455279	228.4354199	252.3437409;
-930.2325581	217.0746395	229.7688152	252.6770119;
-852.7131783	218.4062975	231.1039605	253.0084832;
-775.1937984	219.7389607	232.4413444	253.3419959;
-697.6744186	221.0743056	233.7795541	253.6738734;
-620.1550388	222.4094443	235.1193986	254.0064342;
-542.6356589	223.7437825	236.4626421	254.3388589;
-465.1162791	225.0799554	237.8053721	254.6696124;
-387.5968992	226.4177079	239.150995	255.0008086;
-310.0775194	227.7543394	240.4993028	255.3321426;
-232.5581395	229.0924064	241.8485994	255.6637013;
-155.0387597	230.428923	243.1983212	255.993497;
0	25.66154383	76.54565465	0.039708494;
39.37007874	28.68704857	77.3084298	0.057890648;
78.74015748	31.5326066	78.0722791	0.065700446;
118.1102362	34.27831379	78.8215104	0.063327062;
157.480315	36.91073202	79.57099088	0.052308953;
196.8503937	39.45358777	80.30472701	0.041295248;
236.2204724	41.93994437	81.0305878	0.032008567;
275.5905512	44.34718857	81.74610207	2.47227;
314.9606299	46.7078309	82.45273605	1.98027;
354.3307087	49.02831811	83.14798292	1.7746;
393.7007874	51.30316374	83.83727363	1.92338;
433.0708661	53.55179577	84.52018119	2.51946;
472.4409449	55.76748777	85.19524371	0.036878621;
511.8110236	57.95732438	85.86179966	0.055942822;
551.1811024	60.12396699	86.51955103	0.084546235;
590.5511811	62.26181424	87.18101091	0.12545191;
629.9212598	64.39930194	87.84091725	0.182136641;
669.2913386	66.51663506	88.49723293	0.258904274;
708.6614173	68.62257912	89.16407025	0.360999586;
748.0314961	70.73478117	89.83798005	0.494721689;
787.4015748	72.83397859	90.51727511	0.667532269;
826.7716535	74.94170284	91.21722991	0.88815512;
866.1417323	77.05378811	91.92943758	1.166648602;
905.511811	79.18206124	92.6695248	1.514488774;
944.8818898	81.31316146	93.43060988	1.944713653;
984.2519685	83.46059836	94.22404079	2.472458299;
1023.622047	85.63412741	95.0483985	3.165910552;
1062.992126	87.82141503	95.9137763	3.929208614;
1102.362205	90.02658098	96.81436112	4.846075874;
1141.732283	92.25653028	97.75709486	5.930409822;
1181.102362	94.51178186	98.7418036	7.202702402;
1220.472441	96.78445009	99.77205834	8.6770079;
1259.84252	99.07585298	100.8455801	10.44473275;
1299.212598	101.3887425	101.9628453	12.25928825;
1338.582677	103.711771	103.1177993	14.16376807;
1377.952756	106.0474963	104.312594	16.13303867;
1417.322835	108.3889372	105.5443838	18.1219197;
1456.692913	110.7388413	106.8061514	20.1538726;
1496.062992	113.0849058	108.0918109	22.25905472;
1535.433071	115.4262664	109.4090794	24.37844195;
1574.80315	117.760215	110.741794	26.52619323;
1614.173228	120.088404	112.0907218	28.70085735;
1653.543307	122.3990626	113.450268	30.88619221;
1692.913386	124.6987002	114.8199648	33.1047236;
1732.283465	126.9790018	116.1945894	35.34125003;
1771.653543	129.2473394	117.575273	37.56930983;
1811.023622	131.497214	118.9573292	39.82136358;
1850.393701	133.7279265	120.3395985	42.07621283;
1889.76378	135.9461791	121.7227893	44.32581317;
1929.133858	138.1471551	123.1027072	46.58145653;
1968.503937	140.3390306	124.487629	48.8443321;
2007.874016	142.5158957	125.8761882	51.10206992;
2047.244094	144.6820079	127.2632239	53.36855497;
2086.614173	146.8449869	128.6560895	55.63185074;
2125.984252	148.9999014	130.0589487	57.89636093;
2165.354331	151.1529185	131.4668964	60.15985056;
2204.724409	153.3062529	132.8870286	62.4143256;
2244.094488	155.4615858	134.317878	64.68346774;
2283.464567	157.6165426	135.7649714	66.94856823;
2322.834646	159.7805025	137.2276238	69.22018511;
2362.204724	161.9487033	138.7081442	71.48996501;
2401.574803	164.1243668	140.2098825	73.75869789;
2440.944882	166.3123103	141.7299052	76.03892504;
2480.314961	168.5076614	143.2733879	78.32012817;
2519.685039	170.7124899	144.8348037	80.59724697;
2559.055118	172.9301617	146.4203778	82.88511105;
2598.425197	175.1564645	148.0244647	85.17940827;
2637.795276	177.3971304	149.6516931	87.47056472;
2677.165354	179.6430823	151.2965157	89.77004183;
2716.535433	181.902678	152.9613884	92.06888621;
2755.905512	184.1700737	154.6439782	94.37556656;
2795.275591	186.4445914	156.3410382	96.68323347;
2834.645669	188.7280624	158.0533935	98.99514728;
2874.015748	191.0172945	159.7817527	101.3172611;
2913.385827	193.3132364	161.5219545	103.6380935;
2952.755906	195.6109929	163.2729529	105.9691978;
2992.125984	197.9119781	165.0351838	108.3036169;
3031.496063	200.2154128	166.8108474	110.6515677;
3070.866142	202.5151394	168.5933628	113.000897;
3110.23622	204.8129593	170.38768	115.362273;
3149.606299	207.1027304	172.191967	117.7329594;
3188.976378	209.3807137	174.0024918	120.1190697;
3228.346457	211.6466399	175.8256794	122.5129118;
3267.716535	213.8914991	177.6561261	124.9252503;
3307.086614	216.1115326	179.4980983	127.3501312;
3346.456693	218.3012882	181.348602	129.7888287;
3385.826772	220.4507556	183.2073281	132.2450767;
3425.19685	222.5559467	185.0706441	134.7184362;
3464.566929	224.6058373	186.9438186	137.2027337;
3503.937008	226.5944012	188.8187452	139.7039223;
3543.307087	228.5119106	190.6972009	142.2116103;
3582.677165	230.350164	192.5752722	144.7300099;
3622.047244	232.1016004	194.449507	147.2551269;
3661.417323	233.7600108	196.3188705	149.7813622;
3700.787402	235.3189088	198.1782349	152.3023364;
3740.15748	236.7762304	200.0269772	154.8190637;
3779.527559	238.1299249	201.8594533	157.3228932;
3818.897638	239.3769148	203.6739647	159.8171187;
3858.267717	240.5214455	205.4699692	162.2932413;
3897.637795	241.5645765	207.246589	164.7484701;
3937.007874	242.5120697	209.0016158	167.1853805;
3976.377953	243.3691859	210.7346122	169.6050901;
4015.748031	244.1436521	212.4489239	172.0009738;
4055.11811	244.8422058	214.1431834	174.3778689;
4094.488189	245.4739165	215.8204764	176.739185;
4133.858268	246.0463221	217.4832611	179.0842561;
4173.228346	246.5677309	219.1319776	181.4181974;
4212.598425	247.044939	220.771164	183.7394074;
4251.968504	247.486974	222.4010458	186.0516511;
4291.338583	247.8988242	224.0233829	188.3586793;
4330.708661	248.2850724	225.6430343	190.6635662;
4370.07874	248.6525369	227.2590248	192.9670045;
4409.448819	249.0052353	228.8758622	195.2701406;
4448.818898	249.3439522	230.4922867	197.5772807;
4488.188976	249.6740925	232.1108501	199.8862376;
4527.559055	249.9975336	233.7315811	202.2006129;
4566.929134	250.3154868	235.3549416	204.5206167;
4606.299213	250.6286958	236.9834513	206.8458959;
4645.669291	250.9388681	238.6146352	209.1777373;
4685.03937	251.2453162	240.2517019	211.5165196;
4724.409449	251.5490142	241.8930114	213.8620248;
4763.779528	251.8489536	243.5376694	216.214525;
4803.149606	252.1467172	245.1873052	218.5722081;
4842.519685	252.4397219	246.8397505	220.9393861;
4881.889764	252.7293404	248.4975951	223.3108502;
4921.259843	253.0125665	250.1592147	225.6882838;
4960.629921	253.2925	251.8237762	228.0710638;
5000	253.5647334	253.4919916	230.4600827];
num_lines = 256;
% Interpolate additional values
% Topo_color = interp1(colors_topography.colors_topography(:, 1), colors_topography.colors_topography(:, 2:end), linspace(min(colors_topography.colors_topography(:, 1)), max(colors_topography.colors_topography(:, 1)), num_lines));
Topo_color = interp1(Topography_colors(:, 1),Topography_colors(:, 2:end), linspace(min(Topography_colors(:, 1)), max(Topography_colors(:, 1)), num_lines));

        figure(1);
        set(gcf, 'Units', 'Normalized', 'OuterPosition', [0, 0.04, 1, 0.96],'color','w');clf;
        view(0, 90);
        axesm mercator
        framem;
        gridm;
        ax0 = gca;
        setm(ax0, 'MapProjection', 'robinson');
%         geoshow(Yeq, Xeq,Topo_corrected);
        geoimg = geoshow(Yeq, Xeq,Topo_corrected,'DisplayType', 'texturemap'); % save object
        shadem(-14.5,[210 75]) 
        geoimg.AlphaDataMapping = 'none'; % interpet alpha values as transparency values
        geoimg.FaceAlpha = 'texturemap'; % Indicate that the transparency can be different each pixel
        alpha(geoimg,double(~isnan(Topo_corrected)))
%         ax0.SortMethod = 'childorder';
        
          %Crameri color
%         cmap0=crameri('oleron');
%         colormap(ax0,cmap0);
%         clim([-4000 4000]);

          %ttcmap visualisation
%         [cmap,zlimits] = ttcmap(Topo_corrected,'cmap','mby');
%         colormap(ax0,cmap);
% clim([-8000 4000]);
        
% Number of lines in the gradient
colormap(ax0,Topo_color./num_lines); clim([-8000 4000]);


        c0=colorbar('Position',[0.1 0.18 0.02 0.66]);
        c0.Label.String= "Elevations [m]";
        c0.Label.FontSize = 14;
  
        hold on;
        ax2 = gca;
        setm(ax2, 'MapProjection', 'robinson')
         geoshow(Yeq, Xeq,continents_reversed_position_for_plot_dark,'facealpha', 0);
         contourm(Yeq,Xeq, continents, 0.3, 'k','LineWidth',2);

        contourm(Yeq,Xeq, continents_init, 0.3, 'w','LineWidth',2);
%         quivermc(Yeq_lr, Xeq_lr, V_theta_vec_surface, V_phi_vec_surface,'color','k','linewidth',1,'units','Velocity (m.yr^-^1)','reference',scaling_velocity_factor_surface);
%             plume_str = ['Plumes = ' num2str(numPlumes)];
%             subduction_str = ['Subductions = ' num2str(numSubductions)];
%             title([{time_str} plume_str subduction_str]);
            title({time_str});
        ax2.SortMethod='childorder';
       
        
        %%
        if strcmp(make_video, 'true')
            % Capture the current frame
            frame = getframe(gcf);

            % Write the frame to the video file
            writeVideo(videoWriter, frame);
        else 
            tilefigs;
        end

    end

    if strcmp(make_video, 'true')
        % Close the video file
        close(videoWriter);
    end        
   
end