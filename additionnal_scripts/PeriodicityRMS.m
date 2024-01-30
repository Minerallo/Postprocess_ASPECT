clear all; close all; 

path_model ='/Users/ponsm/Desktop/Kendall_models/v2p3prenew_00_CopyFrederic_114_Disl_1700K_refnu2_100malith_freelslipfreeslip_burnman_phases_pert_coh40_minnu21to300ma_coh20p07_fricdependsonmob_min20_to1ba_rgh_722b4pp1/output/';

depths_for_average = [5000, 15000, 25000, 35000, 45000, 55000, 65000, 75000, 85000, 95000]; %meters
averaged_parameters = 'sinking_velocity';
[time_avg, param_avg] = get_depth_averages(path_model, depths_for_average, averaged_parameters);
% this is velocity, *100 so to be in cm/yr
% y = param_avg{:} .* 100;
for i = 1:numel(param_avg)
    param_avg{i} = param_avg{i} .* 100;
end

% Plotting parameters
line_color = {'r', 'g', 'b', 'k', 'm', 'c', 'y', [0.5 0.5 0.5], [0.5 0.5 0]};
line_style = {'-', '--', '-.', ':'};
marker = {'o', '+', '*', '.', 'x', 's', 'd', '^', 'v', '>', '<', 'p', 'h'};
fig_width = 800;
fig_height = 500;

% Loop through each parameter and plot against time
figure('Position', [100, 100, fig_width, fig_height]);
hold on;
for i = 1:numel(param_avg)
    plot(time_avg, param_avg{i}, ...
         'Color', line_color{mod(i-1, numel(line_color))+1}, ...
         'LineStyle', line_style{mod(i-1, numel(line_style))+1}, ...
         'Marker', marker{mod(i-1, numel(marker))+1},...
         'LineWidth', 1.5);
end
hold off;

% Set axes labels and title
xlabel('Time (years)');
ylabel('Sinking velocity (cm/yr)');
title('Sinking velocity as a function of time');
legend(cellfun(@num2str, num2cell(depths_for_average./1000), 'UniformOutput', false));

% plot periodicity for each param_avg
figure;
hold on;
cmap = jet(numel(param_avg)); % use parula colormap with the same number of colors as the number of depths
for i = 1:numel(param_avg)
    % get data for current parameter
    xxx = time_avg;
    yyy = param_avg{i};

    % find periodicity in data after 0.4e9 yrs
    indexes = find(xxx > 400);
%     if isempty(indexes)
%         continue % skip to next parameter if no data after 0.4e9 yrs
%     end
    minimum_time_row = indexes(1);
    maximum_time_row = indexes(end);

    % only take second half of data due to collapse of lid
    xxx = xxx(minimum_time_row:maximum_time_row);
    yyy = yyy(minimum_time_row:maximum_time_row);

    % calculate the periodicity of the data
    Fs = 1 / (xxx(2) - xxx(1)); % sampling frequency
    L = length(yyy);
    Y = fft(yyy);
    P2 = abs(Y / L);
    P1 = P2(1 : L / 2 + 1);
    P1(2:end-1) = 2 * P1(2:end-1);
    f = Fs * (0 : (L/2)) / L;
    periods = 1 ./ f;
    
    f = f(2:end);
    periods = periods(2:end);

    % truncate periods and Y to have the same length
    n = min(length(periods), length(Y));
    periods = periods(1:n);
    Y = Y(1:n);

    % plot
    scatter(periods, abs(Y)', 100, 'o', 'filled', 'MarkerEdgeColor', 'k', 'MarkerFaceColor', cmap(i,:));
    xlabel('Period (years)');
    ylabel('Power (x10^3)');
    title(['Periodicity vs Time (Period = ', num2str(periods(2)), ' years)']);
    xlim([0, 600]);
%     ylim([0, 0.08]);
end
legend(cellfun(@num2str, num2cell(depths_for_average./1000), 'UniformOutput', false));
