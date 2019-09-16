function cls_out = check_regression_sett(cls_in, def_sett)

% cls_out = check_regression_sett(cls_in, def_sett)
% 
% This function check member variables of chosen classification algorithm 
% structure. It checks if this variables are set correctly or are not set
% and returns the structure with corrected (default) settings.
%
% cls_in            - input classification algorithm structure
% cls_out           - output classification algorithm structure
% def_sett          - default setting for classification algorithms
%
% Supported regression algorithms and settings:
%   a) Classification and Regression Trees algorithm
%       cls.splitcriterion = Split criterion for ClassificationTree
%           {'mse'} 
%           (default: 'mse')
%       cls.minleaf = Minimum number of leaf node observations
%           {1 ... positive integer value}
%           (default: 1
%       cls.minparent = Minimum number of branch node observations
%           {10 ... positive integer value}
%           (default: 10)
%       cls.prune = Pruning flag
%           {'on' 'off'}
%           (default: 'on')
%       cls.prunecriterion = Pruning criterion
%           {'mse'}
%           (default: 'mse)
%
%   b) Ordinal Regression algorithm
%       cls.model = Type of model to fit 
%           {'nominal' 'ordinal' 'hierarchical'}
%           (default: 'ordinal')
%       cls.link  = Link function 
%           {'logit' 'probit' 'comploglog' 'loglog'}
%           (default: 'logit')
%
%   c) Linear Regression algorithm
%       cls.modelspec = Model specification
%           {'constant', 'linear', 'interactions', 'purequadratic', ...
%            'quadratic', 'polyijk'} 
%           (default: 'linear')
%       cls.intercept = Indicator for constant term
%           {'true', 'false'}
%           (default: 'true')
%       cls.robustopts = Indicator of robust fitting type
%           {'on', 'off', string}
%            string: {'andrews', 'bisquare', 'cauchy', 'fair',  ...
%                     'huber', 'logistic', 'ols', 'talwar', 'welsch'}
%           (default (robustopts): 'off')
%           (default (string): 'bisquare')
%
%   d) SVM Regression algorithm
%       cls.kernel = Kernel function
%           {'linear', 'gaussian', 'rbf', 'polynomial'} 
%           (default: 'linear')
%
%   e) Gaussian Process Regression algorithm
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
if ((nargin < 2) || (isempty(def_sett)))  
    def_sett.def_alg                 = 'cart';
    def_sett.cart.alg                = 'cart';
    def_sett.cart.pos_splitcriterion = 'mse';
    def_sett.cart.def_splitcriterion = 'mse';
    def_sett.cart.pos_minleaf        = 1:1:10;
    def_sett.cart.def_minleaf        = 1;
    def_sett.cart.pos_minparent      = 10:1:100;
    def_sett.cart.def_minparent      = 10;
    def_sett.cart.pos_prune          = {'on' 'off'};
    def_sett.cart.def_prune          = 'on';
    def_sett.cart.pos_prunecriterion = 'mse';
    def_sett.cart.def_prunecriterion = 'mse';
    def_sett.or.alg                  = 'or';
    def_sett.or.pos_link             = {'logit' 'probit' 'comploglog' ...
                                        'loglog'};
    def_sett.or.def_link             = 'logit';
    def_sett.or.pos_model            = {'nominal' 'ordinal'           ...
                                        'hierarchical'};
    def_sett.or.def_model            = 'ordinal';
    def_sett.linr.alg                = 'linr';
    def_sett.linr.pos_modelspec      = {'linear', 'constant',         ...
                                        'polyijk' 'interactions',     ...
                                        'purequadratic', 'quadratic'};
    def_sett.linr.def_modelspec      = 'linear';
    def_sett.linr.def_intercept      = {'true' 'false'};
    def_sett.linr.def_intercept      = 'true';
    def_sett.linr.pos_robustopts     = {'on' 'off'};
    def_sett.linr.def_robustopts     = 'off';
    def_sett.svmr.alg                = 'svmr';
    def_sett.svmr.pos_kernel         = {'linear', 'gaussian', 'rbf',  ...
                                        'polynomial'};
    def_sett.svmr.def_kernel         = 'linear';
    def_sett.gpr.alg                 = 'gpr';
    def_sett.gpr.pos_basis           = {'constant', 'none', 'linear', ...
                                        'pureQuadratic'};
    def_sett.gpr.def_basis           = 'constant';
    
else
    
    % Defaults: regression algorithm to use
    if (~isfield(def_sett, 'def_alg'))
        def_sett.def_alg = 'cart';
    end
    
    % Defaults: Classification and Regression Trees settings
    if (~isfield(def_sett, 'cart'))
        def_sett.cart.alg                = 'cart';
        def_sett.cart.pos_splitcriterion = 'mse';
        def_sett.cart.def_splitcriterion = 'mse';
        def_sett.cart.pos_minleaf        = 1:1:10;
        def_sett.cart.def_minleaf        = 1;
        def_sett.cart.pos_minparent      = 10:1:100;
        def_sett.cart.def_minparent      = 10;
        def_sett.cart.pos_prune          = {'on' 'off'};
        def_sett.cart.def_prune          = 'on';
        def_sett.cart.pos_prunecriterion = 'mse';
        def_sett.cart.def_prunecriterion = 'mse';
        
    else
        if (~isfield(def_sett.cart, 'alg'))
            def_sett.cart.alg = 'cart';
        end
        if (~isfield(def_sett.cart, 'pos_splitcriterion'))
            def_sett.cart.pos_splitcriterion = 'mse';
        end
        if (~isfield(def_sett.cart, 'def_splitcriterion'))
            def_sett.cart.def_splitcriterion = 'mse';
        end
        if (~isfield(def_sett.cart, 'pos_minleaf'))
            def_sett.cart.pos_minleaf = 1:1:10;
        end
        if (~isfield(def_sett.cart, 'def_minleaf'))
            def_sett.cart.def_minleaf = 1;
        end
        if (~isfield(def_sett.cart, 'pos_minparent'))
            def_sett.cart.pos_minparent = 10:1:100;
        end
        if (~isfield(def_sett.cart, 'def_minparent'))
            def_sett.cart.def_minparent = 10;
        end
        if (~isfield(def_sett.cart, 'pos_prune'))
            def_sett.cart.pos_prune = {'on' 'off'};
        end
        if (~isfield(def_sett.cart, 'def_prune'))
            def_sett.cart.def_prune = 'on';
        end
        if (~isfield(def_sett.cart, 'pos_prunecriterion'))
            def_sett.cart.pos_prunecriterion = 'mse';
        end
        if (~isfield(def_sett.cart, 'pos_prunecriterion'))
            def_sett.cart.pos_prunecriterion = 'mse';
        end
    end
    
    % Defaults: Ordinal Regression settings
    if (~isfield(def_sett, 'or'))
        def_sett.or.alg       = 'or';
        def_sett.or.pos_link  = {'logit' 'probit' 'comploglog' 'loglog'};
        def_sett.or.def_link  = 'logit';
        def_sett.or.pos_model = {'nominal' 'ordinal' 'hierarchical'};
        def_sett.or.def_model = 'ordinal';
        
    else
        if (~isfield(def_sett.or, 'alg'))
            def_sett.or.alg = 'or';
        end
        if (~isfield(def_sett.or, 'pos_link'))
            def_sett.or.pos_link = {'logit' 'probit' 'comploglog'     ...
                                    'loglog'};
        end
        if (~isfield(def_sett.or, 'def_link'))
            def_sett.or.def_link = 'logit';
        end
        if (~isfield(def_sett.or, 'pos_model'))
            def_sett.or.pos_model = {'nominal' 'ordinal' 'hierarchical'};
        end
        if (~isfield(def_sett.or, 'def_model'))
            def_sett.or.def_model= 'ordinal';
        end
    end
    
    % Defaults: Linear Regression settings
    if (~isfield(def_sett, 'linr'))
        def_sett.linr.alg            = 'linr';
        def_sett.linr.pos_modelspec  = {'linear', 'constant',         ...
                                       'polyijk' 'interactions',      ...
                                       'purequadratic', 'quadratic'};
        def_sett.linr.def_modelspec  = 'linear';
        def_sett.linr.pos_intercept  = {'true', 'false'};
        def_sett.linr.def_intercept  = 'true';
        def_sett.linr.pos_robustopts = {'on', 'off'};
        def_sett.linr.def_robustopts = 'off';
        
    else
        if (~isfield(def_sett.linr, 'alg'))
            def_sett.linr.alg = 'linr';
        end
        if (~isfield(def_sett.linr, 'pos_modelspec'))
            def_sett.linr.pos_modelspec = {'linear', 'constant',      ...
                                           'polyijk' 'interactions',  ...
                                           'purequadratic', 'quadratic'};
        end
        if (~isfield(def_sett.linr, 'def_modelspec'))
            def_sett.linr.def_modelspec = 'linear';
        end
        if (~isfield(def_sett.linr, 'pos_intercept'))
            def_sett.linr.pos_intercept = {'true', 'false'};
        end
        if (~isfield(def_sett.linr, 'def_intercept'))
            def_sett.linr.def_intercept = 'true';
        end
        if (~isfield(def_sett.linr, 'pos_robustopts'))
            def_sett.linr.pos_robustopts = {'on', 'off'};
        end
        if (~isfield(def_sett.linr, 'def_robustopts'))
            def_sett.linr.def_robustopts = 'off';
        end
    end
    
    % Defaults: SVM Regression settings
    if (~isfield(def_sett, 'svmr'))
        def_sett.svmr.alg        = 'svmr';
        def_sett.svmr.pos_kernel = {'linear', 'gaussian', 'rbf',      ...
                                    'polynomial'};
        def_sett.svmr.def_kernel = 'linear';
        
    else
        if (~isfield(def_sett.svmr, 'alg'))
            def_sett.svmr.alg = 'svmr';
        end
        if (~isfield(def_sett.svmr, 'pos_kernel'))
            def_sett.svmr.pos_kernel = {'linear', 'gaussian', 'rbf',  ...
                                        'polynomial'};
        end
        if (~isfield(def_sett.svmr, 'def_kernel'))
            def_sett.svmr.def_kernel = 'linear';
        end
    end
    
    % Defaults: Gaussian Process Regression settings
    if (~isfield(def_sett, 'gpr'))
        def_sett.gpr.alg       = 'gpr';
        def_sett.gpr.pos_basis = {'constant', 'none', 'linear',       ...
                                  'pureQuadratic'};
        def_sett.gpr.def_basis = 'constant';
    
    else
        if (~isfield(def_sett.gpr, 'alg'))
            def_sett.gpr.alg = 'gpr';
        end
        if (~isfield(def_sett.gpr, 'pos_basis'))
            def_sett.gpr.pos_basis = {'constant', 'none', 'linear',   ...
                                      'pureQuadratic'};
        end
        if (~isfield(def_sett.gpr, 'def_basis'))
            def_sett.gpr.def_basis = 'constant';
        end
    end 
end

if (~isfield(cls_in, 'alg'))
    cls_in.alg = def_sett.def_alg;
end

%% Lower the algorithm specification (for consistency)
cls_in.alg = lower(cls_in.alg);

%% Check the classification algorithm description
if (strcmp(cls_in.alg, 'classificationandregressiontrees')         || ...
    strcmp(cls_in.alg, 'classification and regression trees')      || ...
    strcmp(cls_in.alg, 'regressiontrees')                          || ...
    strcmp(cls_in.alg, 'regression trees')                         || ...
    strcmp(cls_in.alg, 'cart algorithm')                           || ...
    strcmp(cls_in.alg, 'cart'))
    
    % Classification and Regression Trees algorithm settings
    def     = def_sett.cart;
    correct = true;
    options = '';
    
    % Check the member variales
    if (~isfield(cls_in, 'splitcriterion'))
        cls_in.splitcriterion = def.def_splitcriterion;
        options = '(default)';
        correct = false;
    else
        if (~sum(strcmpi(cls_in.splitcriterion, def.pos_splitcriterion)))
            cls_in.splitcriterion = def.def_splitcriterion;
            options = '(automatic)';
            correct = false;
        end
    end
    
    if (~isfield(cls_in, 'minleaf'))
        cls_in.minleaf = def.def_minleaf;
        options = '(default)';
        correct = false;
    else
        if (~sum(def.pos_minleaf == cls_in.minleaf))
            cls_in.minleaf = def.def_minleaf;
            options = '(automatic)';
            correct = false;
        end
    end
    
    if (~isfield(cls_in, 'minparent'))
        cls_in.minparent = def.def_minparent;
        options = '(default)';
        correct = false;
    else
        if (~sum(def.pos_minparent == cls_in.minparent))
            cls_in.minparent = def.def_minparent;
            options = '(automatic)';
            correct = false;
        end
    end
    
    if (~isfield(cls_in, 'prune'))
        cls_in.prune = def.def_prune;
        options = '(default)';
        correct = false;
    else
        if (~sum(strcmpi(cls_in.prune, def.pos_prune)))
            cls_in.prune = def.def_prune;
            options = '(automatic)';
            correct = false;
        end
    end
    
    if (~isfield(cls_in, 'prunecriterion'))
        cls_in.prunecriterion = def.def_prunecriterion;
        options = '(default)';
        correct = false;
    else
        if (~sum(strcmpi(cls_in.prunecriterion, def.pos_prunecriterion)))
            cls_in.prunecriterion = def.def_prunecriterion;
            options = '(automatic)';
            correct = false;
        end
    end
    
    cls_out = cls_in;
    
    % Print the corrections or default values (if needed)
%     if (~correct)
%         disp(['Cls settings: 1) Alg: ' cls_in.alg                      ...
%               '; 2) split criterion: ' cls_in.splitcriterion           ...
%               '; 3) minleaf: ' num2str(cls_in.minleaf)                 ...
%               '; 4) minparent: ' num2str(cls_in.minparent)             ...
%               '; 5) prune: ' cls_in.prune                              ...
%               '; 6) pruning criterion: ' cls_in.prunecriterion         ...
%               options]);
%     end
    
elseif (strcmp(cls_in.alg, 'ordinalregression')                     || ...
    strcmp(cls_in.alg, 'ordinal regression')                        || ...
    strcmp(cls_in.alg, 'ordinal reg')                               || ...
    strcmp(cls_in.alg, 'ord regression')                            || ...
    strcmp(cls_in.alg, 'or algorithm')                              || ...
    strcmp(cls_in.alg, 'or'))

    % Ordinal regression algorithm settings
    def     = def_sett.or;
    correct = true;
    options = '';
    
    % Check the member variales
    if (~isfield(cls_in, 'model'))
        cls_in.model = def.def_model;
        options = '(default)';
        correct = false;
    else
        if (~sum(strcmpi(cls_in.model, def.pos_model)))
            cls_in.model = def.def_model;
            options = '(automatic)';
            correct = false;
        end
    end
    
    if (~isfield(cls_in, 'link'))
        cls_in.link = def.def_link;
        options = '(default)';
        correct = false;
    else
        if (~sum(strcmpi(cls_in.link, def.pos_link)))
            cls_in.link = def.def_link;
            options = '(automatic)';
            correct = false;
        end
    end
    
    cls_out = cls_in;
    
    % Print the corrections or default values (if needed)
%     if (~correct)
%         disp(['Cls settings: 1) Alg: ' cls_in.alg                      ...
%               '; 2) model: ' cls_in.model                              ...
%               '; 3) link: ' cls_in.link                                ...
%               options]);
%     end

elseif (strcmp(cls_in.alg, 'linearregression')                      || ...
    strcmp(cls_in.alg, 'linear regression')                         || ...
    strcmp(cls_in.alg, 'linear reg')                                || ...
    strcmp(cls_in.alg, 'lin regression')                            || ...
    strcmp(cls_in.alg, 'linr algorithm')                            || ...
    strcmp(cls_in.alg, 'linr'))

    % Linear regression algorithm settings
    def     = def_sett.linr;
    correct = true;
    options = '';
    
    % Check the member variales
    if (~isfield(cls_in, 'modelspec'))
        cls_in.modelspec = def.def_model;
        options = '(default)';
        correct = false;
    else
        if (~sum(strcmpi(cls_in.modelspec, def.pos_model)))
            cls_in.modelspec = def.def_modelspec;
            options = '(automatic)';
            correct = false;
        end
    end
    
    if (~isfield(cls_in, 'intercept'))
        cls_in.intercept = def.def_intercept;
        options = '(default)';
        correct = false;
    else
        if (~sum(strcmpi(cls_in.intercept, def.pos_intercept)))
            cls_in.intercept = def.def_intercept;
            options = '(automatic)';
            correct = false;
        end
    end
    
    if (~isfield(cls_in, 'robustopts'))
        cls_in.robustopts = def.def_robustopts;
        options = '(default)';
        correct = false;
    else
        if (~sum(strcmpi(cls_in.robustopts, def.pos_robustopts)))
            cls_in.robustopts = def.def_robustopts;
            options = '(automatic)';
            correct = false;
        end
    end
    
    cls_out = cls_in;
    
    % Print the corrections or default values (if needed)
%     if (~correct)
%         disp(['Cls settings: 1) Alg: ' cls_in.alg                      ...
%               '; 2) model: ' cls_in.modelspec                          ...
%               '; 3) intercept: ' cls_in.intercept                      ...
%               '; 4) robustopts: ' cls_in.robustopts                    ...
%               options]);
%     end

elseif (strcmp(cls_in.alg, 'supportvectormachinesregression')       || ...
    strcmp(cls_in.alg, 'support vector machines regression')        || ...
    strcmp(cls_in.alg, 'support vector machines reg')               || ...
    strcmp(cls_in.alg, 'svm regression')                            || ...
    strcmp(cls_in.alg, 'svmr algorithm')                            || ...
    strcmp(cls_in.alg, 'svmr'))

    % SVM regression algorithm settings
    def     = def_sett.svmr;
    correct = true;
    options = '';
    
    % Check the member variales
    if (~isfield(cls_in, 'kernel'))
        cls_in.kernel = def.def_kernel;
        options = '(default)';
        correct = false;
    else
        if (~sum(strcmpi(cls_in.kernel, def.pos_kernel)))
            cls_in.kernel = def.def_kernel;
            options = '(automatic)';
            correct = false;
        end
    end
    
    cls_out = cls_in;
    
    % Print the corrections or default values (if needed)
%     if (~correct)
%         disp(['Cls settings: 1) Alg: ' cls_in.alg                      ...
%               '; 2) kernel: ' cls_in.kernel                            ...
%               options]);
%     end
    
elseif (strcmp(cls_in.alg, 'gaussianprocessregression')             || ...
    strcmp(cls_in.alg, 'gaussian process regression')               || ...
    strcmp(cls_in.alg, 'gaussian process reg')                      || ...
    strcmp(cls_in.alg, 'gp regression')                             || ...
    strcmp(cls_in.alg, 'gpr algorithm')                             || ...
    strcmp(cls_in.alg, 'gpr'))

    % Gaussian process regression algorithm settings
    def     = def_sett.gpr;
    correct = true;
    options = '';
    
    % Check the member variales
    if (~isfield(cls_in, 'basis'))
        cls_in.basis = def.def_basis;
        options = '(default)';
        correct = false;
    else
        if (~sum(strcmpi(cls_in.basis, def.pos_basis)))
            cls_in.basis = def.def_basis;
            options = '(automatic)';
            correct = false;
        end
    end
    
    cls_out = cls_in;
    
    % Print the corrections or default values (if needed)
%     if (~correct)
%         disp(['Cls settings: 1) Alg: ' cls_in.alg                      ...
%               '; 2) basis: ' cls_in.basis                              ...
%               options]);
%     end
    
else
    cls_out = [];
end
