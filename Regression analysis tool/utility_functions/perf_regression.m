function [c_train, c_test, obj] = perf_regression(d_train, ...
    l_train, d_test, cls, mdl)

% [c_train, c_test, obj] = perf_regression(d_train, ...
%     l_train, d_test, cls, mdl)
% 
% The functions performs regression (the training and testing of selected)
% algorithm. The function returns predicted classes for both training and
% testing set of data. The classifier is determined by optional structure
% 'classifier' that holds the information about the algorithm used.
% 
% Supported algorithms:
%   - Classification and Regression trees (CART)
%   - Ordinal Regression (OR)
%   - Linear Regression (LINR)
%   - Support Vector Machines Regression (SVMR)
%   - Gaussian Process Regression (GPR)
%
% INPUT DATA:
% d_train   - input training matrix [columns: features; rows: observations]
% d_test    - input testing matrix  [columns: features; rows: observations]
%             NOTE: size(d_train, 2) must be equal to size(d_test, 2)
% l_train   - column vector with labels of the training observations 
%             NOTE: size(l_train, 1) must be equal to size(d_train, 1)
% cls       - structure with the setting for selected classifier
%             NOTE: .alg field is mandatory (determines the classifier)
%                   possible classifiers: 
%                   - Classification and Regression Trees ('cart')
%                   - Ordinal Regression ('or')
%                   - Linear Regression ('linr')
%                   - Support Vector Machines Regression ('svmr')
%                   - Gaussian Process Regression ('gpr')
% mdl       - trained model object (optional): used in the case when the
%             training process is skipped and the function only performs
%             the regression (testing step)
%
% OUTPUT DATA:
% c_train   - the classified training class (predicted labels)
% c_test    - the classified testing class (predicted labels)
% obj       - trained model object
%
% CLASSIFIER OPTIONAL DATA:
%   1) Classification and Regression trees
%       sett.splitcriterion = Split criterion
%            {'mse'} 
%            (default: 'mse')
%       sett.minleaf = Minimum number of leaf node observations
%            {1 ... positive integer value}
%            (default: 1)
%       sett.minparent = Minimum number of branch node observations
%            {10 ... positive integer value}
%            (default: 10)
%       sett.prune = Pruning flag
%            {'on' 'off'}
%            (default: 'on')
%       sett.prunecriterion = Pruning criterion
%            {'mse'}
%            (default: 'mse)
%
%   2) Ordinal Regression
%       sett.model = Type of model to fit 
%            {'nominal' 'ordinal' 'hierarchical'}
%            (default: 'ordinal')
%       sett.link = Link function 
%            {'logit' 'probit' 'comploglog' 'loglog'}
%            (default: 'logit')
%
%   3) Linear Regression
%       sett.modelspec = Model specification
%            {'constant', 'linear', 'interactions', 'purequadratic', ...
%             'quadratic', 'polyijk'} 
%            (default: 'linear')
%       sett.intercept = Indicator for constant term
%            {'true', 'false'}
%            (default: 'true')
%       sett.robustopts = Indicator of robust fitting type
%            {'on', 'off', string}
%             string: {'andrews', 'bisquare', 'cauchy', 'fair',  ...
%                      'huber', 'logistic', 'ols', 'talwar', 'welsch'}
%            (default (robustopts): 'off')
%            (default (string): 'bisquare')
%
%   4) Support Vector Machine Regression
%       sett.kernel = Kernel function
%            {'linear', 'gaussian', 'rbf', 'polynomial'} 
%            (default: 'linear')
%
%   5) Gaussian Process Regression
%       sett.basis = Explicit basis in the GPR model
%            {'constant', 'none', 'linear', 'pureQuadratic'} 
%            (default: 'constant')
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
if ((nargin < 4) || isempty(cls) || (~isfield(cls, 'alg')))
    error(['The input structure ''cls'' must contain field ''alg'' '   ...
           'that determines the selected regression algorithm.']);
end
if ((nargin < 5) || isempty(mdl))
    mdl = [];
end

%% Lower the algorithm specification (for consistency)
cls.alg = lower(cls.alg);

%% Check the regression algorithm
cls.alg = check_regression_alg(cls.alg);

if (isempty(cls.alg))
    cls.alg = 'cart';
    disp(' ');
    disp('Improper classifier selected. Set to default: Cart');
end

%% Train and classify selected classifier
if (strcmp(cls.alg, 'classificationandregressiontrees')             || ...
    strcmp(cls.alg, 'classification and regression trees')          || ...
    strcmp(cls.alg, 'regressiontrees')                              || ...
    strcmp(cls.alg, 'regression trees')                             || ...
    strcmp(cls.alg, 'cart algorithm')                               || ...
    strcmp(cls.alg, 'cart'))
    
    % Classification and Regression Trees algorithm
    cls = check_regression_sett(cls);
    [c_train, c_test, obj] = ...
        perf_cart(d_train, l_train, d_test, cls, mdl);
    
elseif (strcmp(cls.alg, 'ordinalregression')                        || ...
    strcmp(cls.alg, 'ordinal regression')                           || ...
    strcmp(cls.alg, 'ordinal reg')                                  || ...
    strcmp(cls.alg, 'ord regression')                               || ...
    strcmp(cls.alg, 'or algorithm')                                 || ...
    strcmp(cls.alg, 'or'))
    
    % Ordinal Regression algorithm
    cls = check_regression_sett(cls);
    [c_train, c_test, obj] = ...
        perf_or(d_train, l_train, d_test, cls, mdl);
    
elseif (strcmp(cls.alg, 'linearregression')                         || ...
    strcmp(cls.alg, 'linear regression')                            || ...
    strcmp(cls.alg, 'linear reg')                                   || ...
    strcmp(cls.alg, 'lin regression')                               || ...
    strcmp(cls.alg, 'linr algorithm')                               || ...
    strcmp(cls.alg, 'linr'))
    
    % Linear Regression algorithm
    cls = check_regression_sett(cls);
    [c_train, c_test, obj] = ...
        perf_linr(d_train, l_train, d_test, cls, mdl);

elseif (strcmp(cls.alg, 'supportvectormachinesregression')          || ...
    strcmp(cls.alg, 'support vector machines regression')           || ...
    strcmp(cls.alg, 'support vector machines reg')                  || ...
    strcmp(cls.alg, 'svm regression')                               || ...
    strcmp(cls.alg, 'svmr algorithm')                               || ...
    strcmp(cls.alg, 'svmr'))
    
    % SVM Regression algorithm
    cls = check_regression_sett(cls);
    [c_train, c_test, obj] = ...
        perf_svmr(d_train, l_train, d_test, cls, mdl);
    
elseif (strcmp(cls.alg, 'gaussianprocessregression')                || ...
    strcmp(cls.alg, 'gaussian process regression')                  || ...
    strcmp(cls.alg, 'gaussian process reg')                         || ...
    strcmp(cls.alg, 'gp regression')                                || ...
    strcmp(cls.alg, 'gpr algorithm')                                || ...
    strcmp(cls.alg, 'gpr'))
    
    % Gaussian Process Regression algorithm
    cls = check_regression_sett(cls);
    [c_train, c_test, obj] = ...
        perf_gpr(d_train, l_train, d_test, cls, mdl);
    
else        
    error(['Classifier ' cls.alg ' is not supported.']);
end