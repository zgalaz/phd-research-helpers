function boxPlot(data, options)

% boxPlot(data, options) 
% 
% This function plots the boxplot from input data matrix with the data
% stored in columns. The number of columns is limited to 
% length(options.graph_colors).
%
% -------------------------------------------------------------------------
% INPUT VARIABLES:
% options                    - input settings structure
% data                       - input data values matrix (rows: separate
%                              observations (speakers); cols: data values, 
%                              e.g. speech features 
%
% -------------------------------------------------------------------------
% SETTINGS STRUCTURE:
%   - options.normalize      - data Z-normalization switch [true/false]
%   - options.font_type      - font type to use
%   - options.font_size      - font size to use
%   - options.xlabel         - labels of the x-axis
%   - options.ylabel         - labels of the y-axis
%   - options.labels         - cell array with the labels
%   - options.title          - title of the graph
%   - options.whisker        - maximum whisker length
%
%
%
% --
% ing. Zoltán Galáž
% xgalaz00@stud.feec.vutbr.cz       
% 
% Department of Telecommunications
% Faculty of Electrical Engineering and Communication
% Brno University of Technology

%% Paths and variables
if ((nargin < 2) || isempty(options)) 
    options.normalize = true;
    options.font_type = 'Times New Roman';
    options.font_size = 10;
    options.xlabel    = 'data [-]';
    options.ylabel    = 'data values [-]';
    options.title     = '';
    options.whisker   = 1.5;
    
    num_labels        = 1:size(data, 2);
    options.labels    = cellstr(num2str(num_labels(:)));
else
    if (~isfield(options, 'normalize'))
        options.normalize = true;
    end
    if (~isfield(options, 'font_type'))
        options.font_type = 'Times New Roman';
    end
    if (~isfield(options, 'font_size'))
        options.font_size = 10;  
    end
    if (~isfield(options, 'xlabel'))
        options.xlabel = 'data [-]';
    end
    if (~isfield(options, 'ylabel'))
        options.ylabel = 'data values [-]';
    end
    if (~isfield(options, 'title'))
        options.title = '';
    end
    if (~isfield(options, 'whisker'))
        options.whisker = 1.5;
    end
    if (~isfield(options, 'labels'))
        num_labels     = 1:size(data, 2);
        options.labels = cellstr(num2str(num_labels(:)));
    end
end

%% Set temporary variables (for: code readability)
NUM = size(data, 2);
OPT = options;

%% Normalize feature values
if (OPT.normalize)
    for dat = 1:NUM
        tmp = data(:, dat);
        data(:, dat) = (tmp - mean(tmp))/std(tmp);
    end
end

%% Plot the boxplot
boxplot(data, 'whisker', OPT.whisker, 'labels', OPT.labels); 

%% Set the graph properties
xlabel(OPT.xlabel, 'FontSize', OPT.font_size, 'FontName', OPT.font_type);
ylabel(OPT.ylabel, 'FontSize', OPT.font_size, 'FontName', OPT.font_type); 
title(OPT.title, 'FontSize', OPT.font_size, 'FontName', OPT.font_type);