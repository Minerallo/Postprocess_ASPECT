% Resample the model time step using the resample_topography variable
dt = dt * resample_topography;

% Import the first topography file to determine the model length
filename = fullfile(path_model, 'topography.00008');
topography = import_topography(filename);
modellength = max(topography.x) / 1e3;

% Sort the data in the desired order
x = topography.VarName1;
y = topography.x;
z = topography.y;
r = sqrt(x.^2 + y.^2);
%theta allows to know in which quadrant of the circle we are. 
theta = atan2(y, x); 
%now we can sort the data depending on theta index
[theta_sorted, theta_index] = sort(theta);
x_sorted = x(theta_index);
y_sorted = y(theta_index);
r_sorted = r(theta_index);
z_sorted = z(theta_index);

% Convert theta to degrees and adjust the range to be 0 to 360 degrees
theta_degrees = theta_sorted * 180 / pi;
theta_degrees(theta_degrees<0) = theta_degrees(theta_degrees<0) + 360;

% Find unique theta values and corresponding indices
[theta_unique, index_unique] = unique(theta_degrees);

% Create a new vector for interpolation with high resolution (10000 points)
vec_interp = 0:360/10000:360;

% Interpolate the sorted z data to match the length of theta_degrees
z_interp = interp1(theta_unique, z_sorted(index_unique), vec_interp);



end
