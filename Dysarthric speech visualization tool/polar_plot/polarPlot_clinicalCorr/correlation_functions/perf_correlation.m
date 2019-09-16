function [r, p] = perf_correlation(X, Y, method, options)

% [r, p] = perf_correlation(X, Y, method, options)
% 
% This function computes correlation coefficients and p-values for methods
% specified. It can compute several correlation methods: specified as cell
% array of strings (for several corr. coeff.) or just string (for only one
% corr. coeff. method).
% 
% -----------------------------------------------------------------------
% INPUT VARIABLES:
% X                     - input matrix [observations, instances]
% Y                     - input matrix [observations, instances]
% method                - correlation algorithm specification
% options               - settings (options) of correlation computation
%                        
%                         OPTIONS IS USED ONLY FOR METHODS:
%                         'Spearman', 'Pearson', 'Kendall'
%
%                         MEMBER VARIABLES:
%                           .rows - selected rows for computation
%                                   'all' (default)
%                                   'complete'
%                                   'pairwise'
%                           .tail - the alternative hypothesis against 
%                                   which to compute p-values for testing 
%                                   the hypothesis of no correlation
%                                   'both' (default)
%                                   'right'
%                                   'left'
%
% -----------------------------------------------------------------------
% OUTPUT VARIABLES:
% r                     - correlation coeff.
% p                     - p-values
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
if ((nargin < 4) || (isempty(options)))
    options.rows = 'all';
    options.tail = 'both';
else
    if (~isfield(options, 'rows'))
        options.rows = 'all';
    end
    if (~isfield(options, 'tail'))
        options.tail = 'both';
    end
end

if ((nargin < 3) || (isempty(method)))
    message = ['correlation method uspecified: computing ' ... 
               'Spearman''s corr. coeff., '                ...
               'Pearson''s corr. coeff., '                 ...
               'Kendall''s tau and '                       ...
               'Goodman-Kruskal''s gamma'];
    warning(message);
    method = {'Spearman', 'Pearson', 'Kendal', 'Goodman-Kruskal'};
end

%% Pre-process the methods
if (iscell(method))
    isCellMethod = true;
    numCorrCoeff = length(method);
else
    isCellMethod = false;
    numCorrCoeff = 1;
end

%% Prepare the output variables
r = cell(numCorrCoeff, 1);
p = cell(numCorrCoeff, 1);

%% Compute the corr. coeff.
for coeff = 1:numCorrCoeff
    
    % Set temporary variables
    if (isCellMethod)
        corrCoeff = method{coeff};
    else
        corrCoeff = method;
    end
    
    % Compute the correlations
    if (strcmp(corrCoeff, 'Goodman-Kruskal'))
        corr_r = zeros(size(X, 2), size(Y, 2));
        corr_p = zeros(size(X, 2), size(Y, 2));
        
        for row = 1:size(X, 2)
            for col = 1:size(Y, 2)
                tmp = [X(:, row) Y(:, col)];
                
                [corr_r(row, col), corr_p(row, col)] = gkgammatst(tmp);
            end
        end   
        
    elseif (strcmp(corrCoeff, 'Spearman'))
        [corr_r, corr_p] = corr(X, Y, 'type', 'Spearman', ...
            'rows', options.rows, 'tail', options.tail);
        
    elseif (strcmp(corrCoeff, 'Pearson'))
        [corr_r, corr_p] = corr(X, Y, 'type', 'Pearson',  ...
            'rows', options.rows, 'tail', options.tail);  
        
    elseif (strcmp(corrCoeff, 'Kendall'))
        [corr_r, corr_p] = corr(X, Y, 'type', 'Kendall',  ...
            'rows', options.rows, 'tail', options.tail); 
        
    else
        error('unknown corr. method specified');
    end
    
    % Set the output variables
    r{coeff} = corr_r;
    p{coeff} = corr_p;
end