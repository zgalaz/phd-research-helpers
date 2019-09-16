function kernel_densityPlot(data, options) 

% kernel_densityPlot(data, options)
% 
% This function plots the kernel probability density estimation plot from
% input data matrix with the data stored in columns. The number of columns
% is limited to length(options.graph_colors).
%
% -------------------------------------------------------------------------
% INPUT VARIABLES:
% options                    - input settings structure
% data                       - input data cell array (rows: data values 
%                              (values of speech features); cols: groups
%                              (speaker groups, etc.) 
%
% -------------------------------------------------------------------------
% SETTINGS STRUCTURE:
%   - options.normalize      - data Z-normalization switch [true/false]
%   - options.plot_values    - plot data values switch [true/false]
%   - options.plot_legend    - plot legend switch [true/false]
%   - options.kernel         - kernel function to use
%   - options.npoints        - number of equally spaced points
%   - options.graph_colors   - cell array of graph colors to use
%   - options.graph_lines    - cell array of graph lines to use
%   - options.graph_symbols  - data referencing symbol to use in graphs
%   - options.legend_loc     - legend location to use
%   - options.font_type      - font type to use
%   - options.font_size      - font size to use
%   - options.title          - title of the figure
%   - options.xlabel         - labels of the x-axis
%   - options.ylabel         - labels of the y-axis
%   - options.legend         - cell array with the legend
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
    options.normalize      = true;
    options.plot_values    = true;
    options.plot_legend    = true;
    options.kernel         = 'normal';
    options.npoints        = 100;
    options.graph_colors   = {'b'; 'r'; 'c'; 'm'; 'k'; 'y'};
    options.graph_lines    = {'-'; '--'; '-.'; ':'; '-'; '--'};
    options.graph_symbols  = {'o'; '+'; '*'; 'x'; 's'; 'd'};
    options.legend_loc     = 'northeast';
    options.font_type      = 'Times New Roman';
    options.font_size      = 10;
    options.title          = '';
    options.xlabel         = 'data values [-]';
    options.ylabel         = 'probability density function [-]';
    options.legend         = [];
else
    if (~isfield(options, 'normalize'))
        options.normalize = true;
    end
    if (~isfield(options, 'plot_values'))
        options.plot_values = true;
    end
    if (~isfield(options, 'plot_legend'))
        options.plot_legend = true;
    end
    if (~isfield(options, 'kernel'))
        options.kernel = 'normal';
    end
    if (~isfield(options, 'npoints'))
        options.npoints = 100;
    end
    if (~isfield(options, 'graph_colors'))
        options.graph_colors = {'c'; 'm'; 'b'; 'r'; 'k'; 'y'};
    end
    if (~isfield(options, 'graph_lines'))
        options.graph_lines = {'-'; '--'; '-.'; ':'; '-'; '--'};
    end
    if (~isfield(options, 'graph_symbols'))
        options.graph_symbols = {'o'; '*'; '+'; 'x'; 's'; 'd'};
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
    if (~isfield(options, 'title'))
        options.title = '';  
    end
    if (~isfield(options, 'xlabel'))
        options.xlabel = 'data values [-]';
    end
    if (~isfield(options, 'ylabel'))
        options.ylabel = 'probability density function [-]';
    end
    if (~isfield(options, 'legend'))
        options.legend = [];
    end 
end

%% Check data consistency
if (size(data, 2) > length(options.graph_colors))
    message = ['maximum number of columns exceeded (set to default: ' ...
        num2str(length(options.graph_colors)) ')'];
    warning(message);
    
    % Limit the number of columns in data matrix
    data = data(:, 1:length(options.graph_colors));
end

%% Set temporary variables (for: code readability)
NUM = size(data, 2);
EST = struct();
OPT = options;
brd = ones(1, NUM);
F   = [];
D   = [];
h   = figure;

%% Process the data (per columns)
for dat = 1:NUM
    vec = data(:, dat);
    vec = cell2mat(vec);
    
    % ---------------------------------------------------------------------
    % [A] Normalize feature values
    if (OPT.normalize)
        vec = (vec - mean(vec))/std(vec);
    end
    
    % ---------------------------------------------------------------------
    % [B] Compute the probability density estimation
    [EST(dat).f, EST(dat).x] = ksdensity(vec, ...
        'kernel', OPT.kernel, 'npoints', OPT.npoints);
    
    brd(dat) = length(vec);
    D = [D; vec(:)];
    F = [F EST(dat).f];
end

%% Plot the probability density function estimation
for dat = 1:NUM
    col = [OPT.graph_colors{dat} OPT.graph_lines{dat}];
    plot(EST(dat).x, EST(dat).f, col);
    hold on;
end

%% Plot the actual data values
for dat = 1:NUM
    if (OPT.plot_values)
        
        if (dat == 1)
            low = 1;
        else
            low = brd(dat - 1) + 1;
        end
        
        high = sum(brd(1:dat));
        col  = [OPT.graph_colors{dat} OPT.graph_symbols{dat}];
        
        plot(D(low:high), ones(length(D(low:high)), 1)*(-0.01), col);
        hold on;
    end
end

%% Set the legend
if (OPT.plot_legend)
    if (isempty(OPT.legend))
        L = cell(1, size(data, 2));

        for dat = 1:NUM
            L{dat} = ['data ' num2str(dat)];
        end
    else
        L = OPT.legend;
    end
end

%% Set the graph properties
x_start = min(D(:)) - max(D(:))/2;
y_start = -(max(max(F(:)))/10);
x_end   = max(D(:)) + max(D(:))/2;
y_end   = max(max(F(:))) + max(max(F(:)))/10;

xlim([x_start x_end]); 
ylim([y_start y_end]);

xlabel(OPT.xlabel, 'FontSize', OPT.font_size, 'FontName', OPT.font_type);
ylabel(OPT.ylabel, 'FontSize', OPT.font_size, 'FontName', OPT.font_type);
title(OPT.title, 'FontSize', OPT.font_size, 'FontName', OPT.font_type);

if (OPT.plot_legend)
    legend(L, 'Location', OPT.legend_loc);
end

p = get(h, 'CurrentAxes');
set(p, 'FontSize', OPT.font_size, 'FontName', OPT.font_type);
hold off;
grid off;