function data = load_Damper_data(filename, dataLines)
% If dataLines is not specified, define defaults
if nargin < 2
    dataLines = [2, Inf];
end

%% Set up the Import Options and import the data
opts = delimitedTextImportOptions("NumVariables", 4);

% Specify range and delimiter
opts.DataLines = dataLines;
opts.Delimiter = ";";

% Specify column names and types
opts.VariableNames = ["t (s)"; "p (Pa)"; "dp (Pa)"; "mdot (kg/s)"];
opts.VariableTypes = ["double", "double", "double", "double"];

% Specify file level properties
opts.ExtraColumnsRule = "ignore";
opts.EmptyLineRule = "read";

% Import the data
tmp = readtable(filename, opts);

%% Convert to output type
tmp = table2array(tmp);
data.t = tmp(:,1);
data.p = tmp(:,2);
data.dp = tmp(:,3);
data.mp = tmp(:,4);
end
