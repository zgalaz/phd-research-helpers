function alg = check_regression_alg(label)

% alg = check_regression_alg(label)
% 
% This function check chosen regression algorithm string and returns 
% the string representing chosen algorithm of regression.
%
% label     - regress. algorithm string (one of available possibilities)
% alg       - output string containing the algorithm of regression
%
% Supported algorithms:
%   - Classification and Regression trees (CART)
%   - Ordinal Regression (OR)
%   - Linear Regression (LINR)
%   - Support Vector Machines Regression (SVMR)
%   - Gaussian Process Regression (GPR)
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

%% Lower the algorithm specification (for consistency)
label = lower(label);

%% Check the regression algorithm description
if (strcmp(label, 'classificationandregressiontrees')              || ...
    strcmp(label, 'classification and regression trees')           || ...
    strcmp(label, 'regressiontrees')                               || ...
    strcmp(label, 'regression trees')                              || ...
    strcmp(label, 'cart algorithm')                                || ...
    strcmp(label, 'cart'))
    
    % Classification and Regression Trees algorithm
    alg = 'cart';
    
elseif (strcmp(label, 'ordinalregression')                         || ...
    strcmp(label, 'ordinal regression')                            || ...
    strcmp(label, 'ordinal reg')                                   || ...
    strcmp(label, 'ord regression')                                || ...
    strcmp(label, 'or algorithm')                                  || ...
    strcmp(label, 'or'))
    
    % Ordinal Regression algorithm
    alg = 'or';
    
elseif (strcmp(label, 'linearregression')                          || ...
    strcmp(label, 'linear regression')                             || ...
    strcmp(label, 'linear reg')                                    || ...
    strcmp(label, 'lin regression')                                || ...
    strcmp(label, 'linr algorithm')                                || ...
    strcmp(label, 'linr'))
    
    % Linear Regression algorithm
    alg = 'linr';

elseif (strcmp(label, 'supportvectormachinesregression')           || ...
    strcmp(label, 'support vector machines regression')            || ...
    strcmp(label, 'support vector machines reg')                   || ...
    strcmp(label, 'svm regression')                                || ...
    strcmp(label, 'svmr algorithm')                                || ...
    strcmp(label, 'svmr'))
    
    % SVM Regression algorithm
    alg = 'svmr';
    
elseif (strcmp(label, 'gaussianprocessregression')                 || ...
    strcmp(label, 'gaussian process regression')                   || ...
    strcmp(label, 'gaussian process reg')                          || ...
    strcmp(label, 'gp regression')                                 || ...
    strcmp(label, 'gpr algorithm')                                 || ...
    strcmp(label, 'gpr'))
    
    % Gaussian Process Regression algorithm
    alg = 'gpr';
    
else
    alg = [];
end