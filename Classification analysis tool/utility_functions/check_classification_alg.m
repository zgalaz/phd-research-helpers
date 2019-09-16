function alg = check_classification_alg(label)

% alg = check_classification_alg(label)
% 
% This function check chosen classification algorithm string and returns 
% the string representing chosen algorithm of classification.
%
% label     - classif. algorithm string (one of available possibilities)
% alg       - output string containing the algorithm of classification
%
% Supported algorithms:
%   - Support Vector Machines (SVM)
%   - Naive Bayess (NB)
%   - Linear Discriminat Analysis (LDA)
%   - K-Nearest Neighbor (KNN)
%   - Random Forest (RF)
%   - Gaussian Mixture Models (GMM)
%
%
%
% --
% ing. Zolt�n Gal�
% xgalaz00@stud.feec.vutbr.cz       
% 
% Department of Telecommunications
% Faculty of Electrical Engineering and Communication
% Brno University of Technology

%% Lower the algorithm specification (for consistency)
label = lower(label);

%% Check the classification algorithm description
if (strcmp(label, 'supportvectormachines')                          || ...
    strcmp(label, 'support vector machines')                        || ...
    strcmp(label, 'support vectors')                                || ...
    strcmp(label, 'svmlabel')                                       || ...
    strcmp(label, 'svm label')                                      || ...
    strcmp(label, 'svm'))
    
    % Support Vector Machines label
    alg = 'svm';
        
elseif (strcmp(label, 'naivebayesslabel')                           || ...
        strcmp(label, 'naive bayess label')                         || ...
        strcmp(label, 'naivebayess')                                || ...
        strcmp(label, 'naive bayess')                               || ...
        strcmp(label, 'bayess')                                     || ...
        strcmp(label, 'bayes'))
    
    % Naive Bayess label
    alg = 'bayes';

elseif (strcmp(label, 'lineardiscriminantanalysis')                 || ...
        strcmp(label, 'linear discriminant analysis')               || ...
        strcmp(label, 'lineardiscriminant')                         || ...
        strcmp(label, 'linear discriminant')                        || ...
        strcmp(label, 'lda'))
        
    % Linear Discriminat Analysis label
    alg = 'lda';
    
elseif (strcmp(label, 'knearestneighborlabel')                      || ...
        strcmp(label, 'k-nearestneighborlabel')                     || ...
        strcmp(label, 'k nearest neighbor label')                   || ...
        strcmp(label, 'k-nearest neighbor label')                   || ...
        strcmp(label, 'k nearest neighbor')                         || ...
        strcmp(label, 'k-nearest neighbor')                         || ...
        strcmp(label, 'k-nn')                                       || ...
        strcmp(label, 'knn'))
    
    % K-Nearest Neighbor label
    alg = 'knn';
    
elseif (strcmp(label, 'randomforestslabel')                         || ...
        strcmp(label, 'random forests label')                       || ...
        strcmp(label, 'randomforests')                              || ...
        strcmp(label, 'random forests')                             || ...
        strcmp(label, 'randforests')                                || ...
        strcmp(label, 'rand forests')                               || ...
        strcmp(label, 'rf')                                         || ...
        strcmp(label, 'classificationtreelabel')                    || ...
        strcmp(label, 'classification tree label')                  || ...
        strcmp(label, 'classificationtree')                         || ...
        strcmp(label, 'classification tree')                        || ...
        strcmp(label, 'tree'))
        
    % Random Forest (Classification Tree) label
    alg = 'tree';

elseif (strcmp(label, 'gaussianmixturemodelslabel')                 || ...
        strcmp(label, 'gaussian mixture models label')              || ...
        strcmp(label, 'gaussianmixturemodels')                      || ...
        strcmp(label, 'gaussian mixture models')                    || ...
        strcmp(label, 'gmmlabel')                                   || ...
        strcmp(label, 'gmm label')                                  || ...
        strcmp(label, 'gmm'))
        
    % Gaussian Mixture Models label
    alg = 'gmm';
    
else
    alg = [];
end