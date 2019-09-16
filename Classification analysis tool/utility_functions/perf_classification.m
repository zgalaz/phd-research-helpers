function [c_train, c_test, obj] = perf_classification(d_train, ...
    l_train, d_test, cls, mdl)

% [c_train, c_test, obj] = perf_classification(d_train, ...
%     l_train, d_test, cls, mdl)
% 
% The functions performs classification (training and testing of selected)
% classifier. The function returns predicted classes for both training and
% testing set of data. The classifier is determined by optional structure
% 'classifier' that holds the information about the algorithm used.
% 
% Supported algorithms:
%   - Support Vector Machines (SVM)
%   - Naive Bayess (NB)
%   - Linear Discriminat Analysis (LDA)
%   - K-Nearest Neighbor (KNN)
%   - Random Forest (RF)
%   - Gaussian Mixture Models (GMM)
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
%                   - Support Vector Machines classifier ('svm')
%                   - Naive Bayes classifier ('bayes')
%                   - Linear Discriminant classifier ('lda') 
%                   - K-Nearest Neighbor classifier ('knn') 
%                   - Random Forests classifier ('tree')
%                   - Gaussian Mixture Models classifier ('gmm')
% mdl       - trained model object (optional): used in the case when the
%             training process is skipped and the function only performs
%             the classification (testing step)
%
% OUTPUT DATA:
% c_train   - the classified training class (predicted labels)
% c_test    - the classified testing class (predicted labels)
% obj       - trained model object
% 
% CLASSIFIER OPTIONAL DATA:
%   1) Support Vector Machines classifier
%      cls.kernel = kernel used in SVM algorithm
%                   {'linear' 'quadratic' 'polynomial' 'rbf' 'mlp'}
%                   (default: 'linear')
%
%   2) Naive Bayes classifier
%      cls.distrb = distribution used in Bayes algorithm
%                   {'normal' 'kernel' 'mvmn' 'mn'} 
%                   (default: 'normal')
%      cls.prior  = prior probability for the classes
%                   {'empirical' 'uniform'} 
%                   (default: 'uniform')
% 
%   3) Linear Discriminant classifier
%      cls.prior  = prior probability for the classes
%                   {'empirical' 'uniform'} 
%                   (default: 'uniform')
%      cls.discrm = discrimination type used in the classifier
%                   {'linear' 'quadratic' 'diagLinear' 'diagQuadratic'
%                    'pseudoLinear' 'pseudoQuadratic'} 
%                   (default: 'linear')
% 
%   4) K-Nearest Neighbors classifier
%      cls.k      = the number of nearest neighbors 
%                   (default: 1)
%      cls.dist   = distance calculated in knn classifier
%                   {'euclidean' 'seuclidean' 'cityblock' 'chebychev'
%                    'minkowski' 'mahalanobis' 'cosine' 'correlation'
%                    'spearman' 'hamming' 'jaccard'} 
%                   (default: 'euclidean')
% 
%   5) Random Forests classifier
%      cls.split  = split criterion for ClassificationTree algorithm
%                   {'gdi' 'twoing' 'deviance'} 
%                   (default: 'gdi')
%      cls.prior  = prior probability for the classes
%                   {'empirical', 'uniform'}
%                   (default: 'uniform')
% 
%   6) Gaussian Mixture Models classifier
%      cls.mixt   = number of mixtures of GMM classifier
%                   (default: 1)
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
           'that determines the selected classification algorithm.']);
end
if ((nargin < 5) || isempty(mdl))
    mdl = [];
end

%% Lower the algorithm specification (for consistency)
cls.alg = lower(cls.alg);

%% Check the classification algorithm
cls.alg = check_classification_alg(cls.alg);

if (isempty(cls.alg))
    cls.alg = 'tree';
    disp(' ');
    disp(['Improper classifier selected. Set to default: Tree']);
end

%% Train and classify selected classifier
if (strcmp(cls.alg, 'supportvectormachines')                     || ...
    strcmp(cls.alg, 'support vector machines')                   || ...
    strcmp(cls.alg, 'support vectors')                           || ...
    strcmp(cls.alg, 'svmclassifier')                             || ...
    strcmp(cls.alg, 'svm classifier')                            || ...
    strcmp(cls.alg, 'svm'))
    
    % Support Vector Machines classifier
    cls = check_classification_sett(cls);
    [c_train, c_test, obj] = ...
        perf_svm(d_train, l_train, d_test, cls, mdl);
        
elseif (strcmp(cls.alg, 'naivebayessclassifier')                 || ...
        strcmp(cls.alg, 'naive bayess classifier')               || ...
        strcmp(cls.alg, 'naivebayess')                           || ...
        strcmp(cls.alg, 'naive bayess')                          || ...
        strcmp(cls.alg, 'bayess')                                || ...
        strcmp(cls.alg, 'bayes'))
    
    % Naive Bayess classifier
    cls = check_classification_sett(cls);
    [c_train, c_test, obj] = ...
        perf_bayes(d_train, l_train, d_test, cls, mdl);

elseif (strcmp(cls.alg, 'lineardiscriminantanalysis')            || ...
        strcmp(cls.alg, 'linear discriminant analysis')          || ...
        strcmp(cls.alg, 'lineardiscriminant')                    || ...
        strcmp(cls.alg, 'linear discriminant')                   || ...
        strcmp(cls.alg, 'lda'))
        
    % Linear Discriminat Analysis classifier
    cls = check_classification_sett(cls);
    [c_train, c_test, obj] = ...
        perf_lda(d_train, l_train, d_test, cls, mdl);
    
elseif (strcmp(cls.alg, 'knearestneighborclassifier')            || ...
        strcmp(cls.alg, 'k-nearestneighborclassifier')           || ...
        strcmp(cls.alg, 'k nearest neighbor classifier')         || ...
        strcmp(cls.alg, 'k-nearest neighbor classifier')         || ...
        strcmp(cls.alg, 'k nearest neighbor')                    || ...
        strcmp(cls.alg, 'k-nearest neighbor')                    || ...
        strcmp(cls.alg, 'k-nn')                                  || ...
        strcmp(cls.alg, 'knn'))
    
    % K-Nearest Neighbor classifier
    cls = check_classification_sett(cls);
    [c_train, c_test, obj] = ...
        perf_knn(d_train, l_train, d_test, cls, mdl);
    
elseif (strcmp(cls.alg, 'randomforestsclassifier')               || ...
        strcmp(cls.alg, 'random forests classifier')             || ...
        strcmp(cls.alg, 'randomforests')                         || ...
        strcmp(cls.alg, 'random forests')                        || ...
        strcmp(cls.alg, 'randforests')                           || ...
        strcmp(cls.alg, 'rand forests')                          || ...
        strcmp(cls.alg, 'rf')                                    || ...
        strcmp(cls.alg, 'classificationtreeclassifier')          || ...
        strcmp(cls.alg, 'classification tree classifier')        || ...
        strcmp(cls.alg, 'classificationtree')                    || ...
        strcmp(cls.alg, 'classification tree')                   || ...
        strcmp(cls.alg, 'tree'))
        
    % Random Forest (Classification Tree) classifier
    cls = check_classification_sett(cls);
    [c_train, c_test, obj] = ...
        perf_tree(d_train, l_train, d_test, cls, mdl);

elseif (strcmp(cls.alg, 'gaussianmixturemodelsclassifier')       || ...
        strcmp(cls.alg, 'gaussian mixture models classifier')    || ...
        strcmp(cls.alg, 'gaussianmixturemodels')                 || ...
        strcmp(cls.alg, 'gaussian mixture models')               || ...
        strcmp(cls.alg, 'gmmclassifier')                         || ...
        strcmp(cls.alg, 'gmm classifier')                        || ...
        strcmp(cls.alg, 'gmm'))
        
    % Gaussian Mixture Models classifier
    cls = check_classification_sett(cls);
    [c_train, c_test, obj] = ...
        perf_gmm(d_train, l_train, d_test, cls, mdl);
    
else
    error(['Classifier ' cls.alg ' is not supported.']);
end