function correlationPlot(data, clin, options)

% correlationPlot(data, clin, options)
%
% This function plots so called correlation plot, which is a special type
% of scatter plot in combination with a correlation-based regression line
% showing a correlation relationship between ab input data (e.g. features
% of speech, etc.) and clinical data (e.g. clinical rating scales such as
% UPDRS III score, etc.).
%
% -------------------------------------------------------------------------
% INPUT VARIABLES:
% data                       - input data values (column vector)
% clin                       - input clinical data (column vector)
% options                    - input settings structure
%
% -------------------------------------------------------------------------
% SETTINGS STRUCTURE:
%   - options.type_corr      - type of correlation coeff. to calculate
%   - options.graph_colors   - cell array of graph colors to use
%   - options.graph_symbols  - data referencing symbol to use in graphs
%   - options.num_intervals  - number of intervals to calculate (max: 6)
%   - options.split_function - function to split the data (e.g. @linspace)
%   - options.fit_order      - order of polynomial fit to use
%   - options.legend_loc     - legend location to use
%   - options.font_type      - font type to use
%   - options.font_size      - font size to use
%   - options.xlabel         - labels of the x-axis
%   - options.ylabel         - labels of the y-axis
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
if ((nargin < 3) || isempty(options)) 
    options.type_corr      = 'Spearman';
    options.graph_colors   = {'c'; 'm'; 'b'; 'r'; 'k'; 'y'};
    options.graph_symbols  = '*';
    options.num_intervals  = 5;
    options.split_function = @linspace;
    options.fit_order      = 1;
    options.legend_loc     = 'northeast';
    options.font_type      = 'Times New Roman';
    options.font_size      = 10;
    options.xlabel         = 'independent variable [-]';
    options.ylabel         = 'dependent variable [-]';
else
    if (~isfield(options, 'type_corr'))
        options.type_corr = 'Spearman';  
    end
    if (~isfield(options, 'graph_colors'))
        options.graph_colors = {'c'; 'm'; 'b'; 'r'; 'k'; 'y'};
    end
    if (~isfield(options, 'graph_symbols'))
        options.graph_symbols = '*'; 
    end
    if (~isfield(options, 'num_intervals'))
        options.num_intervals = 5;  
    end
    if (~isfield(options, 'split_function'))
        options.split_function = @linspace;
    end
    if (~isfield(options, 'fit_order'))
        options.fit_order = 1;  
    end
    if (~isfield(options, 'legend_loc'))
        options.legend_loc = 'northeast'; 
    end
    if (~isfield(options, 'font_type'))
        options.font_type = 'Times New Roman';
    end
    if (~isfield(options, 'font_size'))
        options.font_size = 10;  
    end
    if (~isfield(options, 'xlabel'))
        options.xlabel = 'independent variable [-]';
    end
    if (~isfield(options, 'ylabel'))
        options.ylabel = 'dependent variable [-]';
    end
end

%% Check data consistency
if (length(data) ~= length(clin))
    error('observations does not match to clinical data');
end

%% Pre-process the data
delete = ((isnan(data) | isnan(clin)) | ...
          (isinf(data) | isinf(clin)) | ...
          (imag(data)  | imag(clin))) == 1;

data = real(data(~delete));
clin = real(clin(~delete));
OPT  = options;

data = data(:);
clin = clin(:);

%% Calculate the correlation (data / clinical properties)
[r, p] = corr(data, clin, 'type', OPT.type_corr);

%% Prepare the data for the scatter plot
[splits, intervals] = ...
    split_data(clin, OPT.num_intervals, OPT.split_function);

%% Obtain a regression fit (of selected order)
r = round(r*1e4)/1e4;
p = round(p*1e4)/1e4;
q = polyfit(data, clin, OPT.fit_order);

%% Prepare the figure
L = cell(1, OPT.num_intervals);
h = figure;

%% Plot the scatter plot  
for int = 1:OPT.num_intervals

    % Select data from specified interval
    dat = splits(:, int);

    % Set the graph parameters (color, symbol, legend)
    col    = [OPT.graph_symbols OPT.graph_colors{int}];
    L{int} = ['data' ' \in' ' <' num2str(intervals(int)) ...
              ', ' num2str(intervals(int + 1)) ')'];    

    % Plot data from specified interval vs. clinical data 
    plot(data(dat), clin(dat), col);
    hold on;
end

%% Plot the correlation curve
set(refcurve(q), 'Color', 'g', 'LineWidth', 1.5);

%% Set the graph properties
xlabel(OPT.xlabel, 'FontSize', OPT.font_size, 'FontName', OPT.font_type);
ylabel(OPT.ylabel, 'FontSize', OPT.font_size, 'FontName', OPT.font_type);

legend(L, 'Location', OPT.legend_loc);
title(['\rho = ' num2str(r) ', {\itp} = ' num2str(p) ', '   ...
    num2str(OPT.fit_order) '. order fit'],                  ...
    'FontSize', OPT.font_size, 'FontName', OPT.font_type);

p = get(h, 'CurrentAxes');
set(p, 'FontSize', OPT.font_size, 'FontName', OPT.font_type);

axis tight; 
grid on;
hold off;