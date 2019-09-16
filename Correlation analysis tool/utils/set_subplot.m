function [rows, cols, pos] = set_subplot(sett)

% [rows, cols, pos] = set_subplot(sett)
% 
% This function calculates number of rows, number of columns and next 
% position in the subplot object according to input variables, stored
% in the settings structure.
% 
% INPUT STRUCTURE:
% sett.max_columns      - maximum number of columns in subplot 
% sett.num_graphs       - number of columns in subplot to create
% sett.num_features     - number of features to work with (rows)
% sett.actual_pos       - actual position in the subplot object
%
% OUTPUT VARIABLES:
% rows                  - number of rows in the subplot
% cols                  - number of columns in the subplot
% pos                   - next position to use in the subplot

%% Paths and variables
if ((nargin < 1) || (isempty(sett)))
    
    % If the settings structure was not defined, set defaults
    warning('settings for: set_subplot undefined -> use defaults');
    sett.max_columns  = 4;
    sett.num_graphs   = 4;
    sett.num_features = 4;
    sett.actual_pos   = 0;
else
    
    % If the settings structure was not fully defined, set defaults
    if (~isfield(sett, 'max_columns'))
        sett.max_columns = 4;  
        
        str = num2str(sett.max_columns);
        warning(['set_subplot.max_colums undefined -> use ' str]);
    end
    if (~isfield(sett, 'num_graphs'))
        sett.num_graphs = 4;    
        
        str = num2str(sett.num_graphs);
        warning(['set_subplot.num_graphs undefined -> use ' str]);
    end
    if (~isfield(sett, 'num_features'))
        sett.num_features = 4;    
        
        str = num2str(sett.num_features);
        warning(['set_subplot.num_features undefined -> use ' str]);
    end
    if (~isfield(sett, 'actual_pos'))
        sett.actual_pos = 0;  
        
        str = num2str(sett.actual_pos);
        warning(['set_subplot.actual_pos undefined -> use ' str]);
    end
end

%% Set temporary variables (for: readeability)
max_columns  = sett.max_columns;  % maximum number of columns in subplot 
num_graphs   = sett.num_graphs;   % number of columns in subplot to create
num_features = sett.num_features; % number of features to work with (rows)
actual_pos   = sett.actual_pos;   % actual position in the subplot object

%% Set the subplot properties (rows, cols, pos)
if (num_graphs > 1) 
    
    % Description:
    % If number of graphs (e.g. correlation graphs, box plots, etc.) is
    % larger than just one, specified number of rows, number of columns
    % stays the same. The function just increments a position of actual
    % subplot (the output variable "pos" holds the position of the next
    % subplot in the figure)
    
    if (num_graphs > 4)
        error('can not make more than 4 columns in a subplot');
    end
    
    % Inrement the actual subplot position pointer
    if (actual_pos >= 0)
        if (actual_pos < num_features*num_graphs + 1)
            actual_pos = +actual_pos + 1;

            % Set the output variables (rows, cols, pos)
            rows = num_features;
            cols = num_graphs;
            pos  = actual_pos;
        else
            error('max. number of subplot exceeded');
        end
    else
        error('suplots can not contain negative positions')
    end
    
else
     
    % Description:
    % If number of graphs (e.g. correlation graphs, box plots, etc.) is
    % not larger than one, it means we want to plot several features in
    % only one graph representation (e.g. box plots of 12 features), we
    % have to set the number of rows and columns in the subplots to the
    % best fit (max. number of columns is specified in settings struct.)
    
    % Set temporary variable to hold the number of columns needed
    columns = max_columns;
    
    % Decrement the number of columns until it reaches the best fit
    while ((rem(num_features, columns)) ~= 0)
        columns = +columns - 1;
    end
    
    % Inrement the actual subplot position pointer
    if (actual_pos >= 0)
        if (actual_pos < num_features*num_graphs  + 1)
            actual_pos = +actual_pos + 1;

            % Set the output variables (rows, cols, pos)
            rows = num_features/columns;
            cols = columns;
            pos  = actual_pos;
        else
            error('max. number of subplot exceeded');
        end
    else
        error('suplots can not contain negative positions')
    end
end