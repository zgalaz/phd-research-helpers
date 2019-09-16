function [rho, p_val] = correlation_analysis(classification_matrix, ...
                                             feature_labels,        ...
                                             clinical_data,         ...
                                             options) 

% [rho, p_val] = correlation_analysis(classification_matrix, ...
%                                     feature_labels,        ...
%                                     clinical_data,         ...
%                                     options)
% 
% This function calculates a correlation analysis ofor inserted data. It
% computes the chosen correlation coeff. for input classification matrix
% (rows: observations, e.g. patients; columns: features), and resulting
% clinical output, stored in a cell array (rows: observations; columns:
% clinical properties. NOTE: 1. row is determined for the clinical prop.
% names, e.g. UPDRS III, H&Y scale, etc). This function then creates the
% graphs of the features with the best correlation vs. all clinical data
% for all clinical properties stored in "clinical_data" cell array. This
% function also plots the interpolation (regression) line accross these
% plots, using fit with selected order. This function also groups points 
% in the graphs according to severity of disease into chosen number of
% intervals, separated by different colors in the graph. This function 
% returns the matrix of corr. coeffitients (rows: features, columns:
% clinical properties) and p-values matrix, respectively.
%
% INPUTS:
% classification_matrix   - input classification matrix
% feature_labels          - input feature labels (cell array of strings)
% clinical_data           - input clinical data cell array 
% options                 - input settings (options) for the process
%
% OUTPUTS:
% rho                     - output matrix of chosen corr. coeffitients
% p_val                   - output matrix of p-values of correlations
%                           NOTE: rows: features; columns: clin. data
%
% DATA STRUCTURES: 
%   1) classification matrix:
%                                   1. column:       2. column: 
%      1. row: observation n. 1:    feature n. 1:    feature n. 2:                 
%      2. row: observation n. 2:          .                .                
%      .                                  .                .                
%      .                                  .                .                
%
%   2) clinical data:      
%                                   1. column:       2. column:      
%      1. row: clin. data names:    feature n. 1:    feature n. 2:    
%      2. row: clin. data value:          .                .                
%      .                                  .                .                
%      .                                  .                . 
%
% SETTINGS STRUCTURE:
%   - options.type_corr     - type of correlation coeff. to calculate
%   - options.graph_colors  - cell array of graph colors to use
%   - options.graph_symbol  - data referencing symbol to use in graphs
%   - options.num_intervals - number of intervals to calculate (max: 6)
%   - options.fit_order     - order of polynomial fit to use
%   - options.legend_loc    - legend location to use
%   - options.font_type     - font type to use
%   - options.font_size     - font size to use
%   - options.create_xlsx   - switch to create *.xlsx output table

%% Paths and variables
if ((nargin < 4) || isempty(options)) 
    options.type_corr     = 'Spearman';
    options.graph_colors  = {'c'; 'm'; 'b'; 'r'; 'k'; 'y'};
    options.graph_symbol  = '*';
    options.num_intervals = 5;
    options.fit_order     = 3;
    options.legend_loc    = 'northeast';
    options.font_type     = 'Times New Roman';
    options.font_size     = 10;
    options.create_xlsx   = true;
else
    if (~isfield(options, 'type_corr'))
        options.type_corr = 'Spearman';  
    end
    if (~isfield(options, 'graph_colors'))
        options.graph_colors = {'c'; 'm'; 'b'; 'r'; 'k'; 'y'};
    end
    if (~isfield(options, 'graph_symbol'))
        options.graph_symbol = '*'; 
    end
    if (~isfield(options, 'num_intervals'))
        options.num_intervals = 5;  
    end
    if (~isfield(options, 'fit_order'))
        options.fit_order = 3;  
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
    if (~isfield(options, 'create_xlsx'))
        options.create_xlsx = true;  
    end
end

%% Set temporary variables (for: code readability)
CM   = real(classification_matrix);
CLIN = clinical_data;
FL   = feature_labels;
OPT  = options;

%% Obtain input data properties
num_feat = size(CM, 2);   % number of features  (columns)
num_obs  = size(CM, 1);   % number of observations (rows)
num_clin = size(CLIN, 2); % number of clinical features

if ((size(CLIN, 1) - 1) ~= num_obs)
    error('observations does not match to clinical data');
end

%% Prepare the correlation matrices
corr_r = zeros(num_feat, num_clin);
corr_p = zeros(num_feat, num_clin);

%% Calculate the correlation (features / clinical properties)
for clin = 1:num_clin
    
    % Obtain a vector of clinical properties
    vec_clin = real(cell2mat(CLIN(2:end, clin)));

    % Calculate correlation of selected feature with the actual
    % clinical property (feature, e.g. UPDRS III), and store it
    % into outpu matrices (1. rho; 2. p-value)
    for feat = 1:num_feat
        vec_feat = real(CM(:, feat));
        [r, p]   = corr(vec_clin, vec_feat, 'type', OPT.type_corr);
        
        corr_r(feat, clin) = r;
        corr_p(feat, clin) = p;
    end
    
    % Create *.xlsx output file (if selected)
    if (OPT.create_xlsx)
        table = cell(num_feat, length(feature_labels));

        table(:, 1) = FL;
        table(:, 2) = num2cell(corr_r(:, clin));
        table(:, 3) = num2cell(corr_p(:, clin));

        disp(['Storing data to *.xlsx file: ' CLIN{1, clin}]);
        xlswrite('_results\correlation_analysis.xlsx', ...
            table, CLIN{1, clin});
    end
end

%% Plot the correlation graph
h = figure;

num_rows = 3;
num_cols = 3;
for clin = 1:num_clin
    
    % Choose the best matching featue for actual clinical data
    [~, ind] = min(corr_p(:, clin));
    
    % Obtain a vector of clinical properties
    vec_clin = real(cell2mat(CLIN(2:end, clin)));

    % split the observations accorging to clinical properties 
    % into four groups (due to the severity of disorder), e.g. 
    % UPDRS: <0, 25) <25, 50) <50, 75) <75, 100))
    interval_data = false(num_obs, OPT.num_intervals);
    intervals     = linspace(min(vec_clin), ...
                             max(vec_clin), ...
                             OPT.num_intervals + 1);
    
    for int = 1:OPT.num_intervals
        for dat = 1:num_obs
            if ((vec_clin(dat) >= intervals(int)) && ...
                (vec_clin(dat) <  intervals(int + 1)))
                interval_data(dat, int) = true;
            end
        end
    end
    
    % Obtain a vector of observations for the feature with the
    % minimal p-value for chosen correlation coeffitient (rho)
    vec_feat = CM(:, ind);
    
    % Obtain a regression fit (3rd order)
    r = round(corr_r(ind, clin)*1e4)/1e4;
    p = round(corr_p(ind, clin)*1e4)/1e4;
    q = polyfit(vec_feat, vec_clin, OPT.fit_order);
    
    % Obtain a label of the feature with the best correlation
    f_label = FL{ind};
    c_label = CLIN{1, clin};
       
    % -------------------------------------------
    % 1) Scatter plot for feature / clinical data
    L = cell(1, OPT.num_intervals);
    
    for int = 1:OPT.num_intervals
        
        % Select data from specified interval
        dat = interval_data(:, int);
        
        % Set the graph parameters (color, symbol, legend)
        col    = [OPT.graph_symbol OPT.graph_colors{int}];
        L{int} = [c_label ' \in' ' <' num2str(intervals(int)) ...
                  ', ' num2str(intervals(int + 1)) ')'];    
        
        % Plot data from specified interval vs. clinical data
        subplot(num_rows, num_cols, clin);
        plot(vec_feat(dat), vec_clin(dat), col);
        hold on;
    end
    axis tight; 
    grid on;
    
    % Plot reference curve (polynomial fit) and set graph parameters
    rf = refcurve(q);
    set(rf, 'Color', 'g', 'LineWidth', 1.5);
    xlabel(f_label, 'FontSize', OPT.font_size, 'FontName', OPT.font_type);
    ylabel(c_label, 'FontSize', OPT.font_size, 'FontName', OPT.font_type);
    % legend(L, 'Location', OPT.legend_loc);
    title(['\rho = ' num2str(r) ', {\itp} = ' num2str(p)], ...
        'FontSize', OPT.font_size, 'FontName', OPT.font_type);
    
    p = get(h, 'CurrentAxes');
    set(p, 'FontSize', OPT.font_size, 'FontName', OPT.font_type);
    hold off;
end

%% Set the output variables
rho   = corr_r;
p_val = corr_p;
