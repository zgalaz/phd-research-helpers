function score = calc_regression_score(pred, true, alg, data_range, perc)

% score = calc_regression_score(pred, true)
% 
% This function calculates chosen regression score 
%
% pred      - predicted labels (must be numeric)
% true      - true labels (must be numeric)
% alg       - algorithm used to calculate the score (e.g. MAE)
% perc      - return the regression score in percents
%             [0 = OFF, 1 = ON], default: 1 => score [%]
%
% score     - calculated regression score
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

%% Paths and variables / calculate score
if ((nargin < 5) || (isempty(perc)))
    perc = 1;
end
if ((nargin < 4) || (isempty(data_range)))
    data_range = range(true);
end
if ((nargin < 3) || (isempty(alg)))
    alg = 'mae';
end

%% Lower the algorithm string (for consistency)
alg = lower(alg);

%% Calculate selected regression score
if (strcmp(alg, 'ae')                                               || ...
    strcmp(alg, 'abs err')                                          || ...
    strcmp(alg, 'absolute err')                                     || ...
    strcmp(alg, 'absolute error'))

    % Absolute error
    score = absolute_error(true, pred);
    
elseif (strcmp(alg, 'gii')                                          || ...
        strcmp(alg, 'gini')                                         || ...
        strcmp(alg, 'gini idx')                                     || ...
        strcmp(alg, 'gini index'))
    
    % Gini index
    score = gini_index(true, pred);
    
elseif (strcmp(alg, 'mae')                                          || ...
        strcmp(alg, 'mean ae')                                      || ...
        strcmp(alg, 'mean abs err')                                 || ...
        strcmp(alg, 'mean absolute error'))
    
    % Mean absolute error
    score = mean_absolute_error(true, pred);
    
elseif (strcmp(alg, 'mse')                                          || ...
        strcmp(alg, 'mean se')                                      || ...
        strcmp(alg, 'mean sqr err')                                 || ...
        strcmp(alg, 'mean squared error'))
    
    % Mean squared error
    score = mean_squared_error(true, pred);
    
elseif (strcmp(alg, 'msle')                                         || ...
        strcmp(alg, 'mean sle')                                     || ...
        strcmp(alg, 'mean sqr log err')                             || ...
        strcmp(alg, 'mean squared log error'))
    
    % Mean squared log error
    score = mean_squared_log_error(true, pred);
    
elseif (strcmp(alg, 'ngii')                                         || ...
        strcmp(alg, 'ngini')                                        || ...
        strcmp(alg, 'norm gini idx')                                || ...
        strcmp(alg, 'normalized gini index'))
    
    % Normalized Gini index
    score = normalized_gini_index(true, pred);
    
elseif (strcmp(alg, 'rmse')                                         || ...
        strcmp(alg, 'root mean se')                                 || ...
        strcmp(alg, 'root mean sqr err')                            || ...
        strcmp(alg, 'root mean squared error'))
    
    % Root mean squared error
    score = root_mean_squared_error(true, pred);
    
elseif (strcmp(alg, 'rmsle')                                        || ...
        strcmp(alg, 'root mean sle')                                || ...
        strcmp(alg, 'root mean sqr log err')                        || ...
        strcmp(alg, 'root mean squared log error'))
    
    % Root mean squared log error
    score = root_mean_squared_log_error(true, pred);
   
elseif (strcmp(alg, 'se')                                           || ...
        strcmp(alg, 'sqr err')                                      || ...
        strcmp(alg, 'squared err')                                  || ...
        strcmp(alg, 'squared error'))

    % Squared error
    score = squared_error(true, pred);
    
elseif (strcmp(alg, 'sle')                                          || ...
        strcmp(alg, 'sqr log err')                                  || ...
        strcmp(alg, 'squared log err')                              || ...
        strcmp(alg, 'squared log error'))

    % Squared log error
    score = squared_log_error(true, pred);
    
elseif (strcmp(alg, 'eer')                                          || ...
        strcmp(alg, 'rel estim err')                                || ...
        strcmp(alg, 'relative estima err')                          || ...
        strcmp(alg, 'relative estimation error'))

    % Relative estimation error
    score = relative_estimation_error(true, pred, data_range, perc);
    
elseif (strcmp(alg, 'rss')                                          || ...
        strcmp(alg, 'res sum sq')                                   || ...
        strcmp(alg, 'residual sum of sq err')                       || ...
        strcmp(alg, 'residual sum of squares'))

    % Residual sum of squares
    score = residual_sum_of_squares(true, pred);
    
elseif (strcmp(alg, 'cod')                                          || ...
        strcmp(alg, 'coeff of det')                                 || ...
        strcmp(alg, 'coeff of determination')                       || ...
        strcmp(alg, 'coefficient of determination'))

    % Coefficient of determination
    score = coefficient_of_determination(true, pred);
    
else
    error(['Regression score ' alg ' is not supported.']);
end