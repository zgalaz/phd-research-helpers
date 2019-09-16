function cls_out = check_classification_sett(cls_in, def_sett)

% cls_out = check_classification_sett(cls_in, def_sett)
% 
% This function check member variables of chosen classification algorithm 
% structure. It checks if this variables are set correctly or are not set
% and returns the structure with corrected (default) settings.
%
% cls_in            - input classification algorithm structure
% cls_out           - output classification algorithm structure
% def_sett          - default setting for classification algorithms
%
% Supported classification algorithms and settings:
%   a) Support Vector Machines classifier
%       cls.kernel = kernel used in SVM algorithm
%           {'linear' 'quadratic' 'polynomial' 'rbf' 'mlp'}
%           (default: 'linear')
% 
%   b) Naive Bayes classifier
%       cls.distrb = distribution used in Bayes algorithm
%           {'normal' 'kernel' 'mvmn' 'mn'} 
%           (default: 'normal')
%       cls.prior = prior probability for the classes
%           {'empirical' 'uniform'} 
%           (default: 'uniform')
% 
%   c) Linear Discriminant classifier
%       cls.prior = prior probability for the classes
%           {'empirical' 'uniform'} 
%           (default: 'uniform')
%       cls.discrm = discrimination type used in the classifier
%           {'linear' 'quadratic' 'diagLinear' 'diagQuadratic'
%            'pseudoLinear' 'pseudoQuadratic'} 
%           (default: 'linear')
% 
%   d) K-Nearest Neighbors classifier
%       cls.k = the number of nearest neighbors 
%           (default: 1)
%       cls.dist = distance calculated in knn classifier
%           {'euclidean' 'seuclidean' 'cityblock' 'chebychev'
%            'minkowski' 'mahalanobis' 'cosine' 'correlation'
%            'spearman' 'hamming' 'jaccard'} 
%           (default: 'euclidean')
% 
%   e) Random Forests classifier
%       cls.split = split criterion for ClassificationTree alg
%           {'gdi' 'twoing' 'deviance'} 
%           (default: 'gdi')
%       cls.prior = prior probability for the classes
%           {'empirical', 'uniform'}
%           (default: 'uniform')
% 
%   f) Gaussian Mixture Models classifier
%       cls.mixt = number of mixtures of GMM classifier
%           (default: 1)
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
if ((nargin < 2) || (isempty(def_sett)))  
    def_sett.def_alg          = 'tree';
    def_sett.svm.alg          = 'svm';
    def_sett.svm.pos_kernel   = {'linear' 'quadratic' 'polynomial'  ...
                                 'rbf' 'mlp'};
    def_sett.svm.def_kernel   = 'linear';
    def_sett.bayes.alg        = 'bayes';
    def_sett.bayes.pos_distrb = {'normal' 'kernel' 'mvmn' 'mn'};
    def_sett.bayes.def_distrb = 'normal';
    def_sett.bayes.pos_prior  = {'empirical' 'uniform'};
    def_sett.bayes.def_prior  = 'uniform';
    def_sett.lda.alg          = 'lda';
    def_sett.lda.pos_prior    = {'empirical' 'uniform'};
    def_sett.lda.def_prior    = 'uniform';
    def_sett.lda.pos_discrm   = {'linear' 'quadratic' 'diagLinear'  ...
                                 'diagQuadratic' 'pseudoLinear'     ...
                                 'pseudoQuadratic'};
    def_sett.lda.def_discrm   = 'linear';
    def_sett.knn.alg          = 'knn';
    def_sett.knn.pos_k        = 1:10;
    def_sett.knn.def_k        = 1;
    def_sett.knn.pos_dist     = {'euclidean' 'seuclidean'           ...
                                 'cityblock' 'chebychev'            ...
                                 'cosine' 'correlation'             ... 
                                 'spearman' 'minkowski'             ...
                                 'mahalanobis' 'hamming' 'jaccard'};
    def_sett.knn.def_dist     = 'euclidean';
    def_sett.tree.alg         = 'tree';
    def_sett.tree.pos_split   = {'gdi' 'twoing' 'deviance'};
    def_sett.tree.def_split   = 'gdi';
    def_sett.tree.pos_prior   = {'empirical' 'uniform'};
    def_sett.tree.def_prior   = 'uniform';
    def_sett.gmm.alg          = 'gmm';
    def_sett.gmm.pos_mixt     = 1:10;
    def_sett.gmm.def_mixt     = 1;
else
    
    % Default classification algorithm to use
    if (~isfield(def_sett, 'def_alg'))
        def_sett.def_alg = 'tree';
    end
    
    % Default Support Vector Machines algorithm settings
    if (~isfield(def_sett, 'svm'))
        def_sett.svm.alg        = 'svm';
        def_sett.svm.def_kernel = 'linear';
        def_sett.svm.pos_kernel = {'linear' 'quadratic'             ...
                                   'polynomial' 'rbf' 'mlp'};    
    else
        if (~isfield(def_sett.svm, 'alg'))
            def_sett.svm.alg = 'svm';
        end
        if (~isfield(def_sett.svm, 'def_kernel'))
            def_sett.svm.def_kernel = 'linear';
        end
        if (~isfield(def_sett.svm, 'pos_kernel'))
            def_sett.svm.pos_kernel = {'linear' 'quadratic'         ...
                                       'polynomial' 'rbf' 'mlp'};
        end
    end
    
    % Default Naive Bayess algorithm settings
    if (~isfield(def_sett, 'bayes'))
        def_sett.bayes.alg        = 'bayes';
        def_sett.pos_distrb       = {'normal' 'kernel' 'mvmn' 'mn'};
        def_sett.bayes.def_distrb = 'normal';   
        def_sett.bayes.pos_prior  = {'empirical' 'uniform'};
        def_sett.bayes.def_prior  = 'uniform';
    else
        if (~isfield(def_sett.bayes, 'alg'))
            def_sett.bayes.alg = 'bayes';
        end
        if (~isfield(def_sett.bayes, 'pos_distrb'))
            def_sett.bayes.pos_distrb = {'normal' 'kernel' 'mvmn' 'mn'};
        end
        if (~isfield(def_sett.bayes, 'def_distrb'))
            def_sett.bayes.def_distrb = 'normal';
        end
        if (~isfield(def_sett.bayes, 'pos_prior'))
            def_sett.bayes.pos_distrb = {'empirical' 'uniform'};
        end
        if (~isfield(def_sett.bayes, 'def_prior'))
            def_sett.bayes.def_distrb = 'uniform';
        end
    end
    
    % Default Linear Discriminat Analysis algorithm settings
    if (~isfield(def_sett, 'lda'))
        def_sett.lda.alg        = 'lda';
        def_sett.lda.pos_prior  = {'empirical' 'uniform'};
        def_sett.lda.def_prior  = 'uniform';
        def_sett.lda.def_discrm = 'linear';
        def_sett.lda.pos_discrm = {'linear' 'quadratic' 'diagLinear' ...
                                   'diagQuadratic' 'pseudoLinear'    ...
                                   'pseudoQuadratic'};
    else
        if (~isfield(def_sett.lda, 'alg'))
            def_sett.lda.alg = 'lda';
        end
        if (~isfield(def_sett.lda, 'pos_prior'))
            def_sett.lda.pos_prior  = {'empirical' 'uniform'};
        end
        if (~isfield(def_sett.lda, 'def_prior'))
            def_sett.lda.def_prior  = 'uniform';
        end
        if (~isfield(def_sett.lda, 'def_discrm'))
            def_sett.lda.def_discrm = 'linear';
        end
        if (~isfield(def_sett.lda, 'pos_discrm'))
            def_sett.lda.pos_discrm = {'linear' 'quadratic'         ...
                                       'diagLinear' 'diagQuadratic' ...
                                       'pseudoLinear' 'pseudoQuadratic'};
        end
    end
    
    % Default K-Nearest Neighbor algorithm settings
    if (~isfield(def_sett, 'knn'))
        def_sett.knn.alg      = 'knn';
        def_sett.knn.pos_k    = 1:10;
        def_sett.knn.def_k    = 1;
        def_sett.knn.def_dist = 'euclidean';
        def_sett.knn.pos_dist = {'euclidean' 'seuclidean' 'cityblock'  ...
                                 'chebychev' 'minkowski' 'mahalanobis' ...
                                 'cosine' 'correlation' 'spearman'     ...
                                 'hamming' 'jaccard'};
        
    else
        if (~isfield(def_sett.knn, 'alg'))
            def_sett.knn.alg = 'knn';
        end
        if (~isfield(def_sett.knn, 'pos_k'))
            def_sett.knn.pos_k = 1:10;
        end
        if (~isfield(def_sett.knn, 'def_k'))
            def_sett.knn.def_k = 1;
        end
        if (~isfield(def_sett.knn, 'def_dist'))
            def_sett.knn.def_dist = 'euclidean';
        end
        if (~isfield(def_sett.knn, 'pos_dist'))
            def_sett.knn.pos_dist = {'euclidean' 'seuclidean'   ...
                                     'cityblock' 'chebychev'    ...
                                     'minkowski' 'mahalanobis'  ...
                                     'cosine' 'correlation'     ...
                                     'spearman''hamming' 'jaccard'};
        end
    end
    
    % Default Random Forest (Classification Tree) algorithm settings
    if (~isfield(def_sett, 'tree'))
        def_sett.tree.alg       = 'tree';
        def_sett.tree.pos_split = {'gdi' 'twoing' 'deviance'};
        def_sett.tree.def_split = 'gdi';
        def_sett.tree.pos_prior = {'empirical' 'uniform'};
        def_sett.tree.def_prior = 'uniform';
        
    else
        if (~isfield(def_sett.tree, 'alg'))
            def_sett.tree.alg = 'tree';
        end
        if (~isfield(def_sett.tree, 'pos_split'))
            def_sett.tree.pos_split = {'gdi' 'twoing' 'deviance'};
        end
        if (~isfield(def_sett.tree, 'def_split'))
            def_sett.tree.def_split = 'gdi';
        end
        if (~isfield(def_sett.tree, 'pos_prior'))
            def_sett.tree.pos_prior = {'empirical' 'uniform'};
        end
        if (~isfield(def_sett.tree, 'def_prior'))
            def_sett.tree.def_prior = 'uniform';
        end
    end
    
    % Default Gaussian Mixture Models algorithm settings
    if (~isfield(def_sett, 'gmm'))
        def_sett.gmm.alg      = 'gmm';
        def_sett.gmm.pos_mixt = 1:size(tbl_train, 2);
        def_sett.gmm.def_mixt = 1;
        
    else
        if (~isfield(def_sett.gmm, 'alg'))
            def_sett.gmm.alg = 'gmm';
        end
        if (~isfield(def_sett.gmm, 'pos_mixt'))
            def_sett.gmm.pos_mixt = 1:10;
        end
        if (~isfield(def_sett.gmm, 'def_mixt'))
            def_sett.gmm.def_mixt = 1;
        end
    end
end

if (~isfield(cls_in, 'alg'))
    cls_in.alg = def_sett.def_alg;
end

%% Lower the algorithm specification (for consistency)
cls_in.alg = lower(cls_in.alg);

%% Check the classification algorithm description
if (strcmp(cls_in.alg, 'supportvectormachines')                     || ...
    strcmp(cls_in.alg, 'support vector machines')                   || ...
    strcmp(cls_in.alg, 'support vectors')                           || ...
    strcmp(cls_in.alg, 'svmclassifier')                             || ...
    strcmp(cls_in.alg, 'svm classifier')                            || ...
    strcmp(cls_in.alg, 'svm'))
    
    % Support Vector Machines algorithm settings
    def     = def_sett.svm;
    correct = true;
    options = '';
    
    % Check the member variales
    if (~isfield(cls_in, 'kernel'))
        cls_in.kernel = def.def_kernel;
        options       = '(default)';
        correct       = false;
    else
        if (~sum(strcmpi(cls_in.kernel, def.pos_kernel)))
            cls_in.kernel = def.def_kernel;
            options       = '(automatic)';
            correct       = false;
        end
    end
    
    cls_out = cls_in;
    
    % Print the corrections or default values (if needed)
%     if (~correct)
%         disp(['Cls settings: 1) Alg: SVM; 2) kernel: ' ...
%             cls_in.kernel ' ' options]);
%     end
    
elseif (strcmp(cls_in.alg, 'naivebayessclassifier')                 || ...
        strcmp(cls_in.alg, 'naive bayess classifier')               || ...
        strcmp(cls_in.alg, 'naivebayess')                           || ...
        strcmp(cls_in.alg, 'naive bayess')                          || ...
        strcmp(cls_in.alg, 'bayess')                                || ...
        strcmp(cls_in.alg, 'bayes'))
    
    % Naive Bayess algorithm settings
    def       = def_sett.bayes;
    correct   = true;
    options_1 = '';
    options_2 = '';
    
    % Check the member variales
    if (~isfield(cls_in, 'distrb'))
        cls_in.distrb = def.def_distrb;
        options_1     = '(default)';
        correct       = false;
    else
        if (~sum(strcmpi(cls_in.distrb, def.pos_distrb)))
            cls_in.distrb = def.def_distrb;
            options_1     = '(automatic)';
            correct       = false;
        end
    end
    
    if (~isfield(cls_in, 'prior'))
        cls_in.prior = def.def_prior;
        options_2    = '(default)';
        correct      = false;
    else
        if (~sum(strcmpi(cls_in.prior, def.pos_prior)))
            cls_in.prior = def.def_prior;
            options_2    = '(automatic)';
            correct      = false;
        end
    end
    
    cls_out = cls_in;
    
    % Print the corrections or default values (if needed)
%     if (~correct)
%         disp(['Cls settings: 1) Alg: Naive Bayes; 2) distribution: ' ...
%             cls_in.distrb ' ' options_1 '; 3) prior: '               ...
%                 cls_in.prior ' ' options_2]);
%     end

elseif (strcmp(cls_in.alg, 'lineardiscriminantanalysis')            || ...
        strcmp(cls_in.alg, 'linear discriminant analysis')          || ...
        strcmp(cls_in.alg, 'lineardiscriminant')                    || ...
        strcmp(cls_in.alg, 'linear discriminant')                   || ...
        strcmp(cls_in.alg, 'lda'))
        
    % Linear Discriminat Analysis algorithm settings
    def       = def_sett.lda;
    correct   = true;
    options_1 = '';
    options_2 = '';
    
    % Check the member variales
    if (~isfield(cls_in, 'prior'))
        cls_in.prior = def.def_prior;
        options_1    = '(default)';
        correct      = false;
    else
        if (~sum(strcmpi(cls_in.prior, def.pos_prior)))
            cls_in.prior = def.def_prior;
            options_1    = '(automatic)';
            correct      = false;
        end
    end
    
    if (~isfield(cls_in, 'discrm'))
        cls_in.discrm = def.def_discrm;
        options_2     = '(default)';
        correct       = false;
    else
        if (~sum(strcmpi(cls_in.discrm, def.pos_discrm)))
            cls_in.discrm = def.def_discrm;
            options_2     = '(automatic)';
            correct       = false;
        end
    end
    
    cls_out = cls_in;
    
    % Print the corrections or default values (if needed)
%     if (~correct)
%         disp(['Cls settings: 1) Alg: LDA; 2) prior: '           ...
%             cls_in.prior ' ' options_1 '; 3) discrimination: '  ...
%                 cls_in.discrm ' ' options_2]);
%     end
    
elseif (strcmp(cls_in.alg, 'knearestneighborclassifier')            || ...
        strcmp(cls_in.alg, 'k-nearestneighborclassifier')           || ...
        strcmp(cls_in.alg, 'k nearest neighbor classifier')         || ...
        strcmp(cls_in.alg, 'k-nearest neighbor classifier')         || ...
        strcmp(cls_in.alg, 'k nearest neighbor')                    || ...
        strcmp(cls_in.alg, 'k-nearest neighbor')                    || ...
        strcmp(cls_in.alg, 'k-nn')                                  || ...
        strcmp(cls_in.alg, 'knn'))
    
    % K-Nearest Neighbor algorithm settings
    def       = def_sett.knn;
    correct   = true;
    options_1 = '';
    options_2 = '';
    
    % Check the member variales
    if (~isfield(cls_in, 'k'))
        cls_in.k  = def.def_k;
        options_1 = '(default)';
        correct   = false;
    else
        if (~any(cls_in.k == def.pos_k))
            cls_in.k  = def.def_k;
            options_1 = '(automatic)';
            correct   = false;
        end
    end
    
    if (~isfield(cls_in, 'dist'))
        cls_in.dist = def.def_dist;
        options_2   = '(default)';
        correct     = false;
    else
        if (~sum(strcmpi(cls_in.dist, def.pos_dist)))
            cls_in.dist = def.def_dist;
            options_2   = '(automatic)';
            correct     = false;
        end
    end
    
    cls_out = cls_in;
    
    % Print the corrections or default values (if needed)
%     if (~correct)
%         disp(['Cls settings: 1) Alg: K-NN; 2) k: '             ...
%             num2str(cls_in.k) ' ' options_1 '; 3) distance: '  ...
%                 cls_in.dist ' ' options_2]);
%     end

elseif (strcmp(cls_in.alg, 'randomforestsclassifier')               || ...
        strcmp(cls_in.alg, 'random forests classifier')             || ...
        strcmp(cls_in.alg, 'randomforests')                         || ...
        strcmp(cls_in.alg, 'random forests')                        || ...
        strcmp(cls_in.alg, 'randforests')                           || ...
        strcmp(cls_in.alg, 'rand forests')                          || ...
        strcmp(cls_in.alg, 'rf')                                    || ...
        strcmp(cls_in.alg, 'classificationtreeclassifier')          || ...
        strcmp(cls_in.alg, 'classification tree classifier')        || ...
        strcmp(cls_in.alg, 'classificationtree')                    || ...
        strcmp(cls_in.alg, 'classification tree')                   || ...
        strcmp(cls_in.alg, 'tree'))
        
    % Random Forest (Classification Tree) algorithm settings
    def       = def_sett.tree;
    correct   = true;
    options_1 = '';
    options_2 = '';
    
    % Check the member variales
    if (~isfield(cls_in, 'split'))
        cls_in.split = def.def_split;
        options_1    = '(default)';
        correct      = false;
    else
        if (~sum(strcmpi(cls_in.split, def.pos_split)))
            cls_in.split = def.def_split;
            options_1    = '(automatic)';
            correct      = false;
        end
    end
    
    if (~isfield(cls_in, 'prior'))
        cls_in.prior = def.def_prior;
        options_2    = '(default)';
        correct      = false;
    else
        if (~sum(strcmpi(cls_in.prior, def.pos_prior)))
            cls_in.prior = def.def_prior;
            options_2    = '(automatic)';
            correct      = false;
        end
    end
    
    cls_out = cls_in;
    
    % Print the corrections or default values (if needed)
%     if (~correct)
%         disp(['Cls settings: 1) Alg: Tree; 2) split: ' ...
%             cls_in.split ' ' options_1 '; 3) prior: '  ...
%                 cls_in.prior ' ' options_2]);
%     end

elseif (strcmp(cls_in.alg, 'gaussianmixturemodelsclassifier')       || ...
        strcmp(cls_in.alg, 'gaussian mixture models classifier')    || ...
        strcmp(cls_in.alg, 'gaussianmixturemodels')                 || ...
        strcmp(cls_in.alg, 'gaussian mixture models')               || ...
        strcmp(cls_in.alg, 'gmmclassifier')                         || ...
        strcmp(cls_in.alg, 'gmm classifier')                        || ...
        strcmp(cls_in.alg, 'gmm'))
        
    % Gaussian Mixture Models algorithm settings
    def     = def_sett.gmm;
    correct = true;
    options = '';
    
    % Check the member variales
    if (~isfield(cls_in, 'mixt'))
        cls_in.mixt = def.def_mixt;
        options     = '(default)';
        correct     = false;
    else
        if (~any(cls_in.mixt == def.pos_mixt))
            cls_in.mixt = def.def_mixt;
            options     = '(automatic)';
            correct     = false;
        end
    end
    
    cls_out = cls_in;
    
    % Print the corrections or default values (if needed)
%     if (~correct)
%         disp(['Cls settings: 1) Alg: GMM; 2) mixtures: ' ...
%             num2str(cls_in.mixt) ' ' options]);
%     end
    
else
    cls_out = [];
end