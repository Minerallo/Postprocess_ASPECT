function area = polygon_area(lats, lons, radius)
    % polygon_area computes the area of a spherical polygon
    % assuming a spherical Earth.
    % Inputs:
    %   - lats: Latitudes in degrees
    %   - lons: Longitudes in degrees
    %   - algorithm: Unused, included for compatibility
    %   - radius: Radius of the sphere (default is Earth's radius 6378137 meters)
    % Outputs:
    %   - area: Area of the polygon in square meters

    if nargin < 4
        radius = 6378137; % Earth's radius in meters
    end

    % Convert degrees to radians
    lats = deg2rad(lats);
    lons = deg2rad(lons);
    
    % Close the polygon if not already closed
    if lats(1) ~= lats(end)
        lats(end + 1) = lats(1);
        lons(end + 1) = lons(1);
    end
    
    % Calculate colatitudes relative to (0, 0)
    a = sin(lats / 2).^2 + cos(lats) .* sin(lons / 2).^2;
    colat = 2 * atan2(sqrt(a), sqrt(1 - a));
    
    % Calculate azimuths relative to (0, 0)
    az = mod(atan2(cos(lats) .* sin(lons), sin(lats)), 2 * pi);
    
    % Calculate differences in azimuths and colatitudes
    daz = diff(az);
    daz = mod(daz + pi, 2 * pi) - pi;
    deltas = diff(colat) / 2;
    colat = colat(1:end-1) + deltas;
    
    % Perform the line integral to calculate the area
    integrands = (1 - cos(colat)) .* daz;
    
    % Sum the integrals and compute the area
    area_ratio = abs(sum(integrands)) / (4 * pi);
    
    % Ensure the area is the smaller of the two possible values
    area_ratio = min(area_ratio, 1 - area_ratio);
    
    % Return area in units of the radius (default is square meters)
    area = area_ratio * 4 * pi * radius^2;
end
