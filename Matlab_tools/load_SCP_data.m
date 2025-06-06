function data = load_SCP_data(filename, dataLines)
% If dataLines is not specified, define defaults
if nargin < 2
    dataLines = [2, Inf];
end

%% Set up the Import Options and import the data
opts = delimitedTextImportOptions("NumVariables", 10);

% Specify range and delimiter
opts.DataLines = dataLines;
opts.Delimiter = ";";

% Specify column names and types
opts.VariableNames = ["ts", "p0bar", "pLbar", "v0Ms", "vLms", "Q0m3hQLm3hLDLambda", "VarName7", "VarName8", "VarName9", "VarName10"];
opts.VariableTypes = ["double", "double", "double", "double", "double", "double", "double", "double", "double", "double"];

% Specify file level properties
opts.ExtraColumnsRule = "ignore";
opts.EmptyLineRule = "read";

% Import the data
tmp = readtable(filename, opts);

%% Convert to output type
tmp = table2array(tmp);
data.t = tmp(:,1);
data.p_head = tmp(:,2);
data.p_tail = tmp(:,3);
data.v_head = tmp(:,4);
data.v_tail = tmp(:,5);
data.Q_head = tmp(:,6);
data.Q_tail = tmp(:,7);

end
