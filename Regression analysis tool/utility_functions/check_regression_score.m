function score = check_regression_score(label)

% score = check_regression_score(label)
% 
% This function check chosen regression score string and returns 
% the method of regression scoring. 
%
% label     - scoring algorithm string (one of available possibilities)
% score     - output string containing the method of regress. scoring
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

%% Lower the score string (for consistency)
label = lower(label);

%% Set the classification specification (string)
if (strcmp(label, 'ae')                                             || ...
    strcmp(label, 'abs err')                                        || ...
    strcmp(label, 'absolute err')                                   || ...
    strcmp(label, 'absolute error'))

    % Absolute error
    score = 'ae';
    
elseif (strcmp(label, 'gii')                                        || ...
        strcmp(label, 'gini')                                       || ...
        strcmp(label, 'gini idx')                                   || ...
        strcmp(label, 'gini index'))
    
    % Gini index
    score = 'gii';
    
elseif (strcmp(label, 'mae')                                        || ...
        strcmp(label, 'mean ae')                                    || ...
        strcmp(label, 'mean abs err')                               || ...
        strcmp(label, 'mean absolute error'))
    
    % Mean absolute error
    score = 'mae';
    
elseif (strcmp(label, 'mse')                                        || ...
        strcmp(label, 'mean se')                                    || ...
        strcmp(label, 'mean sqr err')                               || ...
        strcmp(label, 'mean squared error'))
    
    % Mean squared error
    score = 'mse';
    
elseif (strcmp(label, 'msle')                                       || ...
        strcmp(label, 'mean sle')                                   || ...
        strcmp(label, 'mean sqr log err')                           || ...
        strcmp(label, 'mean squared log error'))
    
    % Mean squared log error
    score = 'msle';
    
elseif (strcmp(label, 'ngii')                                       || ...
        strcmp(label, 'ngini')                                      || ...
        strcmp(label, 'norm gini idx')                              || ...
        strcmp(label, 'normalized gini index'))
    
    % Normalized Gini index
    score = 'ngii';
    
elseif (strcmp(label, 'rmse')                                       || ...
        strcmp(label, 'root mean se')                               || ...
        strcmp(label, 'root mean sqr err')                          || ...
        strcmp(label, 'root mean squared error'))
    
    % Root mean squared error
    score = 'rmse';
    
elseif (strcmp(label, 'rmsle')                                      || ...
        strcmp(label, 'root mean sle')                              || ...
        strcmp(label, 'root mean sqr log err')                      || ...
        strcmp(label, 'root mean squared log error'))
    
    % Root mean squared log error
    score = 'rmsle';
    
elseif (strcmp(label, 'se')                                         || ...
        strcmp(label, 'sqr err')                                    || ...
        strcmp(label, 'squared err')                                || ...
        strcmp(label, 'squared error'))

    % Squared error
    score = 'se';
    
elseif (strcmp(label, 'sle')                                        || ...
        strcmp(label, 'sqr log err')                                || ...
        strcmp(label, 'squared log err')                            || ...
        strcmp(label, 'squared log error'))

    % Squared log error
    score = 'sle';
    
elseif (strcmp(alg, 'rss')                                          || ...
        strcmp(alg, 'res sum sq')                                   || ...
        strcmp(alg, 'residual sum of sq err')                       || ...
        strcmp(alg, 'residual sum of squares'))
    
    % Residual sum of squares
    score = 'rss';
    
elseif (strcmp(alg, 'cod')                                          || ...
        strcmp(alg, 'coeff of det')                                 || ...
        strcmp(alg, 'coeff of determination')                       || ...
        strcmp(alg, 'coefficient of determination'))
    
    % Coefficient of determination
    score = 'cod';
    
else
    score = [];
end