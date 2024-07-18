clc; 
clear all; 
close all;

% Define the paths to the images for each model
imagePath1 = '/Users/poulami/Nextcloud/3D_spherical_model/model_global_Poulami_postprocess/14_test_250Ma_3D_test_with_shear_heat_01/Maps_geofeatures';
imagePath2 = '/Users/poulami/Nextcloud/3D_spherical_model/model_global_Poulami_postprocess/13_test_250Ma_3D_test_no_shear_heat/Maps_geofeatures';
imagePath3 = '/Users/poulami/Nextcloud/3D_spherical_model/model_global_Poulami_postprocess/12_test_250Ma_3D_test_no_shear_heat/Maps_geofeatures';
imagePath4 = '/Users/poulami/Nextcloud/3D_spherical_model/model_global_Poulami_postprocess/11_test_250Ma_3D_test_no_shear_heat/Maps_geofeatures';
imagePath5 = '/Users/poulami/Nextcloud/3D_spherical_model/model_global_Poulami_postprocess/03_visc_cmb_extended_15_3D_no_shear_heat/Maps_geofeatures';

% Define the depth slice (you need to specify this value)
depthSlice = '2800'; % Example depth slice

model_titles={'title 1';'title 2';'title 3';'title 4';'title 5'};
% Define the times in Myr for all models
times = [1, 51, 101, 151, 201, 251];

% Number of time points
numTimes = length(times);

% Define the cropping rectangle [xmin ymin width height]
cropRect = [500 300 2500 1150];

% Create a figure with white background
figure('Color', 'white');

% Calculate subplot positions
subplotHeight = 0.13; % Height for each subplot
spacing = 0.005; % Spacing between subplots (adjust as needed)
width_subplot = 0.16; % Adjust width of each subplot column

% Main title height
mainTitleHeight = 0.07;

% Loop through each model
for model = 1:5
    % Calculate the position for the main title of this model
    mainTitlePosition = [0.1 + (width_subplot + spacing) * (model - 1), 1 - mainTitleHeight, width_subplot, mainTitleHeight];
    
    % Create a subplot for the main title of this model without frame
    ax_mainTitle = axes('Position', mainTitlePosition, 'Color', 'none', 'XTick', [], 'YTick', [], 'Box', 'off', 'Visible', 'off');
    text(0.5, 0.5, sprintf(model_titles{model}), 'FontSize', 14, ...
        'HorizontalAlignment', 'center', 'VerticalAlignment', 'middle');
    
    % Loop through each time and plot the corresponding images for this model
    for i = 1:numTimes
        % Calculate the position for the subplot of this model and time
        subplotBottom = 1 - (i * (subplotHeight + spacing)) - mainTitleHeight; % Bottom position for each subplot
        subplotPosition = [0.1 + (width_subplot + spacing) * (model - 1), subplotBottom, width_subplot, subplotHeight]; % [left, bottom, width, height]
        
        % Create a subplot for this model and time image
        ax = axes('Position', subplotPosition);
        
        % Remove grid from subplot
        grid(ax, 'off');
        
        % Construct the filename for this model and time
        filename = sprintf('geofeatures_depth_%s_%04d.png', depthSlice, times(i));
        
        % Full path to the image for this model and time
        switch model
            case 1
                imagePathFull = fullfile(imagePath1, filename);
            case 2
                imagePathFull = fullfile(imagePath2, filename);
            case 3
                imagePathFull = fullfile(imagePath3, filename);
            case 4
                imagePathFull = fullfile(imagePath4, filename);
            case 5
                imagePathFull = fullfile(imagePath5, filename);
        end
        
        % Read the image for this model and time
        img = imread(imagePathFull);
        
        % Crop the image for this model and time
        imgCropped = imcrop(img, cropRect);
        
        % Display the cropped image for this model and time
        imshow(imgCropped, 'Parent', ax);
if model==1
    % Add time text on the left side of the subplot for the first model
    text(ax, -0.1, 0.5, sprintf('%d Myr', times(i)), 'FontSize', 14, ...
        'HorizontalAlignment', 'right', 'VerticalAlignment', 'middle', ...
        'Units', 'normalized');
end
    end
end

% Adjust figure properties
set(gcf, 'Position', [0, 0, 1500, 800]); % Adjust size as needed
