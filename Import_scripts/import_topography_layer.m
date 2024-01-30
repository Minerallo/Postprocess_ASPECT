function topolayer = import_topography_layer(filename)
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
            topolayer = readtable(filename, opts);

end