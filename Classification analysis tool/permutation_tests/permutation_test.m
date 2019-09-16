function [p_value, out] = permutation_test(tbl_train, label_train, sett)

% [p_value, out] = permutation_test(tbl_train, label_train, sett)
% 
% This function performs the permutation test: by randomly permuting the 
% observations across classes and calculating classi?cation score (such as
% classification accuracy) at each permutation it is possible to establish
% an empirical null distribution of classi?cation accuracies on the random 
% observations. Therefore, the tails of this distribution can then be used
% to determine signi?cance boundaries for a given rate of tolerated false 
% positives (i.e. correct classi?cations that occur by chance).
%
%       ------------------------------------------------------------------
%       Required inputs:
% 
%           tbl_train:          Input feature matrix. Columns: features; 
%                               rows: observations. The matrix is numeric
%                               and it represent the training and testing
%                               set of data for the permutation process.
%
%           label_train:        Input labels that corresponds to classes
%                               of the observations from tbl_train. This
%                               algorithm requires labels to be numeric,
%                               where: 1 - positive hypothesis (disease)
%                                      0 - negative hypothesis (healthy)
%
%       ------------------------------------------------------------------
%       Optional inputs:
%
%           num_cross_val.:     The number of cross validation runs that 
%                               the algorithm perform in every permutation 
%                               run. The score after the num_cross_val is 
%                               consequently averaged in order to obtain 
%                               estimation of the classification score in
%                               the cross validation fold
%                               Default: 100 cross-validation runs
%
%           k_fold:             The number of cross-validation folds in 
%                               every cross-validation run
%                               Default: 10-fold cross-validation
%
%           classifier:         The classifier to use in oreder to do the
%                               classification of the testing data. This
%                               is specified in a structure holding assoc.
%                               member variables settings for the cls
%                               Default: .alg   - 'tree'
%                                        .split - 'gdi'
%                                        .prior - 'uniform'
%
%           score:              The scoring algorithm that is used for the
%                               permutation testing purposes to asses the
%                               statistical significance level of proposed
%                               classifieation model
%                               Default: 'acc' - classification accuracy
%
%           num_permutations:   The number of permutation runs to perform
%                               Defalur: 10 000
%
%       ------------------------------------------------------------------
%       Outputs:
%
%           p_value:            p value of the classification model after
%                               the permutation test of num_permutations
%
%           out:                Structure holding the output variables:
%                               out.base 
%                               out.perm
%                               out.sett                    Enumeration:
%                               ------------------------------------------
%                               out.sett.num_cross_validations
%                               out.sett.num_permutations
%                               out.sett.k_fold
%                               out.sett.classifier
%                               out.sett.score
%                               out.base.k_fold 
%                               out.base.k_fold.contingency_table
%                               out.base.score
%                               out.perm.k_fold
%                               out.perm.k_fold.contingency_table 
%                               out.perm.score
%                               out.perm.p_value
%
% Implemented according to:
% Combrisson E, Jerbi K. Exceeding chance level by chance: The caveat of
% theoreti-cal chance levels in brain signal classification and statistical
% assessment of decoding accuracy. J Neurosci Methods (2015)
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

%% Define the optional variable settings
%  [A] Define all the possible classification algorithms and also its
%      associated member variables. The supported classifiers are listed
%      bellow with the description of settings + default values selected
%      if the option of the classifier is not specified, or is specified
%      incorrectly (corrections are displayed into the console)

%      Classification algorithms and settings:
%           a) Support Vector Machines classifier
%              cls.kernel = kernel used in SVM algorithm
%                   {'linear' 'quadratic' 'polynomial' 'rbf' 'mlp'}
%                   (default: 'linear')
%
%           b) Naive Bayes classifier
%               cls.distrb = distribution used in Bayes algorithm
%                   {'normal' 'kernel' 'mvmn' 'mn'} 
%                   (default: 'normal')
%               cls.prior  = prior probability for the classes
%                   {'empirical' 'uniform'} 
%                   (default: 'uniform')
% 
%           c) Linear Discriminant classifier
%               cls.prior = prior probability for the classes
%                   {'empirical' 'uniform'} 
%                   (default: 'uniform')
%               cls.discrm = discrimination type used in the classifier
%                   {'linear' 'quadratic' 'diagLinear' 'diagQuadratic'
%                    'pseudoLinear' 'pseudoQuadratic'} 
%                   (default: 'linear')
% 
%           d) K-Nearest Neighbors classifier
%               cls.k = the number of nearest neighbors 
%                   (default: 1)
%               cls.dist = distance calculated in knn classifier
%                    {'euclidean' 'seuclidean' 'cityblock' 'chebychev'
%                     'minkowski' 'mahalanobis' 'cosine' 'correlation'
%                     'spearman' 'hamming' 'jaccard'} 
%                    (default: 'euclidean')
% 
%           e) Random Forests classifier
%               cls.split = split criterion for ClassificationTree alg
%                   {'gdi' 'twoing' 'deviance'} 
%                   (default: 'gdi')
%               cls.prior = prior probability for the classes
%                   {'empirical', 'uniform'}
%                   (default: 'uniform')
% 
%           f) Gaussian Mixture Models classifier
%               cls.mixt = number of mixtures of GMM classifier
%                   (default: 1)

classifier.def_alg          = 'tree';
classifier.pos_alg          = {'tree' 'svm' 'bayes' 'lda' 'knn' 'gmm'};
classifier.svm.alg          = 'svm';
classifier.svm.pos_kernel   = {'linear' 'quadratic' 'polynomial' ...
                               'rbf' 'mlp'};
classifier.svm.def_kernel   = 'linear';
classifier.bayes.alg        = 'bayes';
classifier.bayes.pos_distrb = {'normal' 'kernel' 'mvmn' 'mn'};
classifier.bayes.def_distrb = 'normal';
classifier.bayes.pos_prior  = {'empirical' 'uniform'};
classifier.bayes.def_prior  = 'uniform';
classifier.lda.alg          = 'lda';
classifier.lda.pos_prior    = {'empirical' 'uniform'};
classifier.lda.def_prior    = 'uniform';
classifier.lda.pos_discrm   = {'linear' 'quadratic' 'diagLinear' ...
                               'diagQuadratic' 'pseudoLinear'    ...
                               'pseudoQuadratic'};
classifier.lda.def_discrm   = 'linear';
classifier.knn.alg          = 'knn';
classifier.knn.pos_k        = 1:size(tbl_train, 2);
classifier.knn.def_k        = 1;
classifier.knn.pos_dist     = {'euclidean' 'seuclidean' 'cityblock'  ...
                               'chebychev' 'minkowski' 'mahalanobis' ...
                               'cosine' 'correlation' 'spearman'     ...
                               'hamming' 'jaccard'};
classifier.knn.def_dist     = 'euclidean';
classifier.tree.alg         = 'tree';
classifier.tree.pos_split   = {'gdi' 'twoing' 'deviance'};
classifier.tree.def_split   = 'gdi';
classifier.tree.pos_prior   = {'empirical' 'uniform'};
classifier.tree.def_prior   = 'uniform';
classifier.gmm.alg          = 'gmm';
classifier.gmm.pos_mixt     = 1:size(tbl_train, 2);
classifier.gmm.def_mixt     = 1;

%  [B] Define the classification scoring algorithm supported. The list of
%      all algorithms is listed bellow. Default scoring algorithm is used
%      if input structure does not contain information about algorithm to
%      use or if it is incorrectly inserted. Default algorithm: 'acc'
%
%      Scoring algorithms (used in binary classification):
%           - acc (classification accuracy)            
%           - f1  (classification f1 score)
%           - inf (classification informedness)
%           - mar (classification markedness)
%           - pre (classification prevalence)
%           - sen (classification sensitivity)
%           - spe (classification sensitivity)
%           - dor (diagnostic odds ratio)
%           - fdr (false discovery rate)
%           - fnr (false negative rate)
%           - for (false omission rate)
%           - fpr (false positive rate)
%           - mcc (Mathews correlation coeff.)
%           - nlr (negative likelihood ratio)
%           - npv (negative predictive value)
%           - plr (positive likelihood ratio)
%           - ppv (positive predictive value)
%           - tss (tradeoff between SEN and SPE)
%           - yjs (Youden's J statistics)

score.pos_alg = {{'acc' '(classification accuracy)'};         ...
                 {'f1'  '(classification f1 score)'};         ...
                 {'inf' '(classification informedness)'};     ...
                 {'mar' '(classification markedness)'};       ...
                 {'pre' '(classification prevalence)'};       ...
                 {'sen' '(classification sensitivity)'};      ...
                 {'spe' '(classification specificity)'};      ...
                 {'dor' '(diagnostic odds ratio)'};           ...
                 {'fdr' '(false discovery rate)'};            ...
                 {'fnr' '(false negative rate)'};             ...
                 {'for' '(false omission rate)'};             ...
                 {'fpr' '(false positive rate)'};             ...
                 {'mcc' '(Mathews correlation coeff.)'};      ...
                 {'nlr' '(negative likelihood ratio)'};       ...
                 {'npv' '(negative predictive value)'};       ...
                 {'plr' '(positive likelihood ratio)'};       ...
                 {'ppv' '(positive predictive value)'};       ...
                 {'yjs' '(Youden''s J statistics)'};          ...
                 {'tss' '(tradeoff between SEN and SPE)'}};
score.def_alg =  'acc';

%% Set the default/helper variables
%  [A] Helper variables that allow optional variables check and a console 
%      output. Those variables contain the helper variables (set / unset)
%      and default variables that hold the default values for unspecified
%      variables (empty, not specified or out of range). The algorithm is
%      implemented to apply the predefined value

DEF_NUM_CROSS_VAL    = true; % Helper var.: number of CV rounds specified
DEF_K_FOLD           = true; % Helper var.: k (number) CV folds specified
DEF_CLASSIFIER       = true; % Helper var.: classifier algorithm specified
DEF_SCORE            = true; % Helper var.: cls scoring algorithm specified
DEF_NUM_PERMUTATIONS = true; % Helper var.: num. permutations specified       
DEF_PERMUTATIONS_MTX = true; % Helper var.: mtx. of permutations specified 

%  [B] Default settings for the SFFS algorithm. It includes the setting of
%      the number of cross-validation runs and k-folds of cross-validation
%      then the classifier algorithm, scoring function, number of features
%      to look ahead without closing the algorithm and the maximum number
%      of features to be inside the selected feature set

DEFAULT_NUM_CROSS_VAL    = 100;                % DEFAULT: num. of CV rounds 
DEFAULT_K_FOLD           = 10;                 % DEFAULT: nun. of CV folds
DEFAULT_CLASSIFIER       = classifier.def_alg; % DEFAULT: cls algorithm
DEFAULT_SCORE            = score.def_alg;      % DEFAULT: scoring algorithm
DEFAULT_NUM_PERMUTATIONS = 10000;              % DEFAULT: num. permutations
DEFAULT_DISP_STAT        = false;              % DEFAULT: disp. statistics
DEFAULT_DISP_HIST        = false;              % DEFAULT: disp. histogram
DEFAULT_DISP_P_PROGRESS  = false;              % DEFAULT: disp. histogram
DEFAULT_PERMUTATIONS_MTX = false;              % DEFAULT: mtx. of perms

%% Paths and variables
if ((nargin < 3) || (isempty(sett)))  
    
    % Set the number of cross-validations
    sett.num_cross_validations = DEFAULT_NUM_CROSS_VAL;
    DEF_NUM_CROSS_VAL = false;
    
    % Set the k-folds of cross-validation 
    sett.k_fold = DEFAULT_K_FOLD;
    DEF_K_FOLD = false;
    
    % Set the classification algorithm
    sett.classifier = struct();
    sett.classifier.alg = DEFAULT_CLASSIFIER;
    DEF_CLASSIFIER = false;
    
    % Set the scoring function of classification
    sett.score = DEFAULT_SCORE;
    DEF_SCORE = false;
    
    % Set the number of permutations
    sett.num_permutations = DEFAULT_NUM_PERMUTATIONS;
    DEF_NUM_PERMUTATIONS = false;
    
    % Set the matrix of permutations
    sett.permutations_mtx = DEFAULT_PERMUTATIONS_MTX;
    DEF_PERMUTATIONS_MTX = false;
    
    % Set the switch for the properties displaying 
    sett.disp_stat = DEFAULT_DISP_STAT;
    
    % Set the switch for the histogram displaying 
    sett.disp_hist = DEFAULT_DISP_HIST;
    
    % Set the switch for the p-value progress displaying 
    sett.disp_p_progress = DEFAULT_DISP_P_PROGRESS;
    
else
    
    % Set the number of cross-validations
    if (~isfield(sett, 'num_cross_validations'))
        sett.num_cross_validations = DEFAULT_NUM_CROSS_VAL;
        DEF_NUM_CROSS_VAL = false; 
    end
    
    % Set the k-folds of cross-validation 
    if (~isfield(sett, 'k_fold'))
        sett.k_fold = DEFAULT_K_FOLD;
        DEF_K_FOLD = false;
    end
    
    % Set the classification algorithm 
    if (~isfield(sett, 'classifier'))
        sett.classifier = struct();
        sett.classifier.alg = DEFAULT_CLASSIFIER;
        DEF_CLASSIFIER = false;
    end
    
    % Set the scoring function of classification
    if (~isfield(sett, 'score'))
        sett.score = DEFAULT_SCORE;
        DEF_SCORE = false;
    end
    
    % Set the number of permutations
    if (~isfield(sett, 'num_permutations'))
        sett.num_permutations = DEFAULT_NUM_PERMUTATIONS;
        DEF_NUM_PERMUTATIONS = false;
    end
    % Set the matrix of permutations
    if (~isfield(sett, 'permutations_mtx'))
        sett.permutations_mtx = DEFAULT_PERMUTATIONS_MTX;
        DEF_PERMUTATIONS_MTX = false;
    end
    
    % Set the switch for the properties displaying 
    if (~isfield(sett, 'disp_stat'))
        sett.disp_stat = DEFAULT_DISP_STAT;
    end
    
    % Set the switch for the histogram displaying 
    if (~isfield(sett, 'disp_hist'))
        sett.disp_hist = DEFAULT_DISP_HIST;
    end
    
    % Set the switch for the p-value progress displaying 
    if (~isfield(sett, 'disp_p_progress'))
        sett.disp_p_progress = DEFAULT_DISP_P_PROGRESS;
    end
end

sett.score = lower(sett.score);
sett.classifier.alg = lower(sett.classifier.alg);

%% Test the input data   
%  [A] Check if the number of observations in 'tbl_train' is equal to the
%      number of observations in 'label_train'. The algorithm demands the 
%      same number of observations (rows) = e.g. patients (PD/HC, etc.)

if (size(tbl_train, 1) ~= size(label_train, 1))
    error('The number of observations in tbl_train ~= label_train');
end

%  [B] Optional variables handling routines to prepare the console output
%
%      Handle: The number of cross-validation runs
%              DEFAULT: 100 cross-validation runs
if (~DEF_NUM_CROSS_VAL)
    STR_NUM_CROSS_VAL = num2str(sett.num_cross_validations);
    STR_NUM_CROSS_VAL_STATE = ' (DEFAULT)';
elseif ((DEF_NUM_CROSS_VAL) && (sett.num_cross_validations <= 0))
    sett.num_cross_validations = DEFAULT_NUM_CROSS_VAL;
    STR_NUM_CROSS_VAL = num2str(sett.num_cross_validations);
    STR_NUM_CROSS_VAL_STATE = ' (AUTOMATIC)';
else
    STR_NUM_CROSS_VAL = num2str(sett.num_cross_validations);
    STR_NUM_CROSS_VAL_STATE = '';
end

%      Handle: The number of cross-validation folds (k)
%              DEFAULT: 10 cross-validation folds
if (~DEF_K_FOLD)
    STR_K_FOLD = [num2str(sett.k_fold) '-fold'];
    STR_K_FOLD_STATE = ' (DEFAULT)';
elseif ((DEF_K_FOLD) && ((sett.k_fold <= 0) || ...
        (sett.k_fold > floor(size(tbl_train, 1)/2))))
    sett.k_fold = DEFAULT_K_FOLD;
    STR_K_FOLD = [num2str(sett.k_fold) '-fold'];
    STR_K_FOLD_STATE = ' (AUTOMATIC)';
elseif ((DEF_K_FOLD) && (sett.k_fold == 1))    
    STR_K_FOLD = 'L-O-OUT';
    STR_K_FOLD_STATE = '';

    if (sett.num_cross_validations ~= 1)
        sett.num_cross_validations = 1;
        STR_NUM_CROSS_VAL = num2str(sett.num_cross_validations);
        STR_NUM_CROSS_VAL_STATE = ' (AUTOMATIC)';
    end
else       
    STR_K_FOLD = [num2str(sett.k_fold) '-fold'];
    STR_K_FOLD_STATE = '';
end

%      Handle: Classification algorithm used in the feature selection
%              DEFAULT: ClassificatioTree: (Random Forest) classifier
if (~DEF_CLASSIFIER)
    STR_CLASSIFIER = sett.classifier.alg;
    STR_CLASSIFIER_STATE = ' (DEFAULT)';
    sett.classifier = ...
            check_classification_sett(sett.classifier, classifier);
else

    % check the classifier (set the classification algorithm)
    sett.classifier.alg = check_classification_alg(sett.classifier.alg);
    if (~isempty(sett.classifier.alg))
         STR_CLASSIFIER = sett.classifier.alg;
         STR_CLASSIFIER_STATE = '';
         sett.classifier = ...
            check_classification_sett(sett.classifier, classifier);
    else
        disp(' '                                        );
        disp('Unknown classifier, set to DEFAULT: tree' );
        disp('Known classifiers:'                       );
        disp('   - tree  (Classification Tree)'         );
        disp('   - svm   (Support Vector Machines)'     );
        disp('   - knn   (K-Nearest Neighbors)'         );
        disp('   - gmm   (Gaussian Mixture Models)'     );
        disp('   - lda   (Linear discriminant Analysis)');
        disp('   - bayes (Naive Bayes Classifier)'      );
        
        STR_CLASSIFIER = DEFAULT_CLASSIFIER;
        STR_CLASSIFIER_STATE = ' (AUTOMATIC)';
        sett.classifier = ...
            check_classification_sett(sett.classifier, classifier);
    end
end

%      Handle: The ccoring method used used in the feature selection
%              DEFAULT: 'tss' (tradeoff between SEN and SPE)
if (~DEF_SCORE)
    STR_SCORE = sett.score;
    STR_SCORE_STATE = ' (DEFAULT)';
else

    % check the classification score (set the classification score)
    sett.score = check_classification_score(sett.score);
    if (~isempty(sett.score))
         STR_SCORE = sett.score;
         STR_SCORE_STATE = '';
    else
        disp(' '                                           );
        disp('Unknown scoring method, set to DEFAULT: tss' );
        disp('Known scoring methods:'                      );
        disp('   - acc (classification accuracy)'          );            
        disp('   - f1  (classification f1 score)'          );
        disp('   - inf (classification informedness)'      );
        disp('   - mar (classification markedness)'        );
        disp('   - pre (classification prevalence)'        );
        disp('   - sen (classification sensitivity)'       );
        disp('   - spe (classification sensitivity)'       );
        disp('   - dor (diagnostic odds ratio)'            );
        disp('   - fdr (false discovery rate)'             );
        disp('   - fnr (false negative rate)'              );
        disp('   - for (false omission rate)'              );
        disp('   - fpr (false positive rate)'              );
        disp('   - mcc (Mathews correlation coeff.)'       );
        disp('   - nlr (negative likelihood ratio)'        );
        disp('   - npv (negative predictive value)'        );
        disp('   - plr (positive likelihood ratio)'        );
        disp('   - ppv (positive predictive value)'        );
        disp('   - tss (Youden''s J statistics)'           );
        disp('   - yjs (tradeoff between SEN and SPE)'     );
        
        STR_SCORE = DEFAULT_SCORE;
        STR_SCORE_STATE = ' (AUTOMATIC)';
    end
end

%      Handle: Number of permutations used to assess the classifier
%              DEFAULT: 10000
if (~DEF_NUM_PERMUTATIONS)
    STR_NUM_PERMUTATIONS = num2str(sett.num_permutations);
    STR_NUM_PERMUTATIONS_STATE = ' (DEFAULT)';
else
    if (sett.num_permutations <= 0)
        str = '10000';
        disp(' '                                                 );
        disp(['Improper num. permutations, set to DEFAULT: ' str]);

        sett.num_permutations = DEFAULT_NUM_PERMUTATIONS;
        STR_NUM_PERMUTATIONS = num2str(sett.num_permutations);
        STR_NUM_PERMUTATIONS_STATE = ' (AUTOMATIC)';
    else
        STR_NUM_PERMUTATIONS = num2str(sett.num_permutations);
        STR_NUM_PERMUTATIONS_STATE = '';
    end
end

%      Handle: Matrix of permutations used to assess the classifier
%              DEFAULT: perm per iter
if (~DEF_PERMUTATIONS_MTX)
    STR_PERMUTATIONS_MTX = 'perm per iter';
    STR_PERMUTATIONS_MTX_STATE = ' (DEFAULT)';
else
    STR_PERMUTATIONS_MTX = 'defined';
    STR_PERMUTATIONS_MTX_STATE = '';
end

%  [C] Console output of the input variables and warnings (if occurred)
if (sett.disp_stat)
    fprintf('\n'                                                    );
    fprintf('INPUT VARIABLES\n'                                     );
    fprintf('training matrix:\t\t(OK): specified\n'                 );
    fprintf('training labels:\t\t(OK): specified\n'                 );
    fprintf('\n'                                                    );
    fprintf('OPTIONAL VARIABLES\n'                                  );
    fprintf('Cross validations:\t\t(OK): %s%s\n',  ...
        STR_NUM_CROSS_VAL, STR_NUM_CROSS_VAL_STATE                  );
    fprintf('K-fold       \t\t\t(OK): %s%s\n',     ...
        STR_K_FOLD, STR_K_FOLD_STATE                                );
    fprintf('Classifier:\t\t\t\t(OK): %s%s\n',     ...
        STR_CLASSIFIER, STR_CLASSIFIER_STATE                        );
    fprintf('Scoring method:\t\t\t(OK): %s%s\n',   ...
        STR_SCORE, STR_SCORE_STATE                                  );
    fprintf('Num. permutations:\t\t(OK): %s%s\n',   ...
        STR_NUM_PERMUTATIONS, STR_NUM_PERMUTATIONS_STATE            );
    fprintf('Mtx. permutations:\t\t(OK): %s%s\n',   ...
        STR_PERMUTATIONS_MTX, STR_PERMUTATIONS_MTX_STATE            );
end

%% Set the possible classes (check classes)
if ((iscell(label_train)) || (length(unique(label_train)) > 2))
    error('Labels must be a numerical vector with < 3 classes');
else
    ncv = sett.num_cross_validations;
    cls = sett.classifier;
    
    if (sett.k_fold == 1)
        fld = size(tbl_train, 1);
    else
        fld = sett.k_fold;
    end
    
    out      = struct();
    out.base = struct();
    out.perm = struct();
    out.sett = struct();
    out.base.num_cv.cv_fold.contingency_table = zeros(2, 2);
    out.perm.num_cv.cv_fold.contingency_table = zeros(2, 2);
    out.perm.p_value = 0;
    
    out.sett.num_cross_validations = sett.num_cross_validations;
    out.sett.num_permutations = sett.num_permutations;
    out.sett.k_fold = sett.k_fold;
    out.sett.score = sett.score;
    out.sett.classifier = sett.classifier;
    
    scr_perm    = zeros(sett.num_permutations, 1);    
    pred_labels = cell(ncv, fld);
    true_labels = cell(ncv, fld);
end

%% Calculate the baseline classification score
%
%  Perform the cross validations (number of validations is specified by 
%  the settings of the algorithm). In each of the validation steps data 
%  is randomly shuffled and next set (training/testing) data's selected 
%  (stratification)
for cvi = 1:ncv
                
    % Prepare cross validation data and labels for the
    % cross validation. From this set the k-folds will 
    % be executed. It also randomizes the data
    if (sett.k_fold == 1)
        cv = cvpartition(numel(label_train), 'LeaveOut');
    else
        cv = cvpartition(numel(label_train), ...
            'KFold', sett.k_fold);
    end

    % Create temporary vector (predictions/presennted)
    predictions = cell(1, fld);
    presented   = cell(1, fld);

    % Iterate over all folds of the cross validation
    % and do the training and testing process of the
    % chosen rclassification algorithm
    for fold = 1:fld

        % Set the indices for the cross-validation
        train_idx = cv.training(fold);
        test_idx  = cv.test(fold);

        % Perform the estimation   
        train_tbl = tbl_train(train_idx, :);
        train_lbl = label_train(train_idx);
        test_tbl  = tbl_train(test_idx, :);
        test_lbl  = label_train(test_idx);

        [~, pred_nested, cls_model] = ...
            perf_classification(train_tbl, train_lbl, ...
                test_tbl, cls);

        predictions{fold} = pred_nested(:);
        presented{fold}   = test_lbl(:);
    end

    % Update the true and predicted labels
    pred_labels(cvi, :) = predictions;
    true_labels(cvi, :) = presented;
end

%% Compute overall classification performance score that is achieved
SCR_cv  = zeros(ncv, fld);
for cvi = 1:ncv
    for fold = 1:fld

        true_D = (true_labels{cvi, fold} == 1);
        detc_D = (pred_labels{cvi, fold} == 1);
        true_H = (true_labels{cvi, fold} == 0);
        detc_H = (pred_labels{cvi, fold} == 0);

        % Calculate the values of the confusion matrix
        TP = sum(true_D & detc_D);
        TN = sum(true_H & detc_H);
        FP = sum(true_H & detc_D);
        FN = sum(true_D & detc_H);

        tbl = [TP FP; FN TN];
        out.base.num_cv(cvi).cv_fold(fold).contingency_table = tbl;
        
        % Calculate the classification scores
        SCR_cv(cvi, fold) = calc_classifiation_score( ...
            TP, FP, FN, TN, sett.score);
    end 
end

out.base.score = mean(SCR_cv(:));

%% Perform the permutation test for presented classification model 
for p = 1:sett.num_permutations

    % Check if the matrix of permutations is/is not defined
    if (~DEF_PERMUTATIONS_MTX)
        rng(p, 'twister'); 
        perm_labels = randperm(length(label_train)).';      
    else
        perm_labels = sett.permutations_mtx(:, p);
    end
     
    % Permute the labels before the cross validation
    label_train_p = label_train(perm_labels, :);
    
    % Perform the cross validation process
    for cvi = 1:ncv
        
        % Prepare cross validation data and labels for the
        % cross validation. From this set the k-folds will 
        % be executed. It also randomizes the data
        if (sett.k_fold == 1)
            cv = cvpartition(numel(label_train_p), 'LeaveOut');
        else
            cv = cvpartition(numel(label_train_p), ...
                'KFold', sett.k_fold);
        end

        % Create temporary vector (predictions/presennted)
        predictions = cell(1, fld);
        presented   = cell(1, fld);

        % Iterate over all folds of the cross validation
        % and do the training and testing process of the
        % chosen rclassification algorithm
        for fold = 1:fld

            % Set the indices for the cross-validation
            train_idx = cv.training(fold);
            test_idx  = cv.test(fold);

            % Perform the estimation   
            train_tbl = tbl_train(train_idx, :);
            train_lbl = label_train_p(train_idx);
            test_tbl  = tbl_train(test_idx, :);
            test_lbl  = label_train_p(test_idx);

            [~, pred_nested, cls_model] = ...
                perf_classification(train_tbl, train_lbl, ...
                    test_tbl, cls);

            predictions{fold} = pred_nested(:);
            presented{fold}   = test_lbl(:);
        end

        % Update the true and predicted labels
        pred_labels(cvi, :) = predictions;
        true_labels(cvi, :) = presented;
    end
    
    % Compute overall classification performance score that is achieved
    % so far (default: 'acc' class. scoring alg.) 
    SCR_cv  = zeros(ncv, fld);
    for cvi = 1:ncv
        for fold = 1:fld

            true_D = (true_labels{cvi, fold} == 1);
            detc_D = (pred_labels{cvi, fold} == 1);
            true_H = (true_labels{cvi, fold} == 0);
            detc_H = (pred_labels{cvi, fold} == 0);

            % Calculate the values of the confusion matrix
            TP = sum(true_D & detc_D);
            TN = sum(true_H & detc_H);
            FP = sum(true_H & detc_D);
            FN = sum(true_D & detc_H);

            tbl = [TP FP; FN TN];
            out.perm(p).num_cv(cvi).cv_fold(fold).contingency_table = tbl;
            
            % Calculate the classification scores
            SCR_cv(cvi, fold) = calc_classifiation_score( ...
                TP, FP, FN, TN, sett.score);
        end 
    end

    % Calculate the mean scores over all cross validations
    scr_perm(p) = mean(SCR_cv(:));
    
    out.perm(p).score   = scr_perm(p);
    out.perm(p).p_value = sum(out.base.score < scr_perm(1:p))/p;
end

%% Set the output variables
p_value = sum(out.base.score < scr_perm)/sett.num_permutations;

%% Plot the histogram (if specified)
if (sett.disp_hist)
    figure;
    
    b = 0:10/sett.num_permutations:1;
    h = histc(scr_perm, b);
    bar(b, h);
    hold on;
    line([out.base.score out.base.score], [0 max(h)]);
    hold off;
    title(sprintf('scr = %.3f', out.base.score));
end

%% Plot the p-value progress (if specified)
if (sett.disp_p_progress)
    perm_vec = zeros(1, sett.num_permutations);
    
    for p = 1:sett.num_permutations
        perm_vec(p) = out.perm(p).p_value;
    end
    
    figure;
    plot(1:sett.num_permutations, perm_vec, 'r');
    hold on; grid on;
    line([1 sett.num_permutations], [p_value p_value]);
    hold off;
    title(sprintf('p-value progress: %d iterations', p));
    xlabel('\rightarrow {\it number of permutations} [-]');
    ylabel('\rightarrow {\it p-value} [-]');
end