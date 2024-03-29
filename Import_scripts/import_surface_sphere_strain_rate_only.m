function surface = import_surface_sphere_strain_rate_only(filename, dataLines)
%IMPORTFILE Import data from a text file
%  SURFACE0000 = IMPORTFILE(FILENAME) reads data from text file FILENAME
%  for the default selection.  Returns the data as a table.
%
%  SURFACE0000 = IMPORTFILE(FILE, DATALINES) reads data for the
%  specified row interval(s) of text file FILENAME. Specify DATALINES as
%  a positive scalar integer or a N-by-2 array of positive scalar
%  integers for dis-contiguous row intervals.
%
%  Example:
%  surface0000 = importfile("/Users/ponsm/Desktop/Model_sphere/Models_HLRN/R01a_Rodinia_1GPa_Mantle_C40MPa_LR/surface_0000.csv", [2, Inf]);
%
%  See also READTABLE.
%
% Auto-generated by MATLAB on 24-Jan-2024 00:23:47

%% Input handling

% If dataLines is not specified, define defaults
if nargin < 2
    dataLines = [2, Inf];
end

%% Set up the Import Options and import the data
opts = delimitedTextImportOptions("NumVariables", 10);

% Specify range and delimiter
opts.DataLines = dataLines;
opts.Delimiter = ",";

% Specify column names and types
opts.VariableNames = ["Time", "Points0", "Points1", "Points2", "continent", "strain_rate", "v_phi", "v_r", "v_theta", "vtkValidPointMask"];
opts.VariableTypes = ["double", "double", "double", "double", "double", "double", "double", "double", "double", "double"];

% Specify file level properties
opts.ExtraColumnsRule = "ignore";
opts.EmptyLineRule = "read";

% Import the data
surface = readtable(filename, opts);

end