function data = load_LWP_data(filename, dataLines)
% If dataLines is not specified, define defaults
if nargin < 2
    dataLines = [2, Inf];
end

%% Set up the Import Options and import the data
opts = delimitedTextImportOptions("NumVariables", 11);

% Specify range and delimiter
opts.DataLines = dataLines;
opts.Delimiter = ";";

% Specify column names and types
opts.VariableNames = ["t (s)", "p(0) (bar)", "p(L) (bar)","v(0) m/s", "v(L) (m/s)","T(0) (C)"," T(L) (C)","rho(0) (kg/m3)"," rho(L) (kg/m3)","mp(0) (kg/s)"," mp(L) (kg/s)"];
opts.VariableTypes = ["double", "double", "double", "double", "double", "double", "double", "double", "double", "double", "double"];

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
data.T_head = tmp(:,6);
data.T_tail = tmp(:,7);
data.rho_head = tmp(:,8);
data.rho_tail = tmp(:,9);
data.mp_head = tmp(:,10);
data.mp_tail = tmp(:,11);
end
