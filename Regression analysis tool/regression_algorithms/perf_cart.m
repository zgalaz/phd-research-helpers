function [c_train, c_test, obj] = perf_cart(d_train, l_train, ...
    d_test, sett, mdl)

% [c_train, c_test, obj] = perf_cart(d_train, l_train, d_test, sett)
% 
% This function performs regression (training; testing) of Classification
% and Regression Trees algorithm. This function returns predicted classes
% (multi-classes => prediction of 'continuous' values) for both training 
% and testing set of data. The classifier is  specified by setings (sett)
% structure.
% 
% INPUT DATA:
% d_train   - input training matrix [columns: features; rows: observations]
% d_test    - input testing matrix  [columns: features; rows: observations]
%             NOTE: size(d_train, 2) must be equal to size(d_test, 2)
% l_train   - column vector with labels of the training observations 
%             NOTE: size(l_train, 1) must be equal to size(d_train, 1)
% sett      - structure with the setting for selected classifier
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
% for more information, see: 
%   <http://www.mathworks.com/help/stats/fitrtree.html>
%   <http://www.mathworks.com/help/stats/regressiontree-class.html>
%
% sett.splitcriterion = Split criterion for ClassificationTree algorithm
%      {'mse'} 
%      (default: 'mse')
% sett.minleaf = Minimum number of leaf node observations
%      {1 ... positive integer value}
%      (default: 1)
% sett.minparent = Minimum number of branch node observations
%      {10 ... positive integer value}
%      (default: 10)
% sett.prune = Pruning flag
%      {'on' 'off'}
%      (default: 'on')
% sett.prunecriterion = Pruning criterion
%      {'mse'}
%      (default: 'mse)
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
if ((nargin < 4) || (isempty(sett)))
    sett.splitcriterion  = 'mse';
    sett.minleaf         = 1;
    sett.minparent       = 10;
    sett.prune           = 'on';
    sett.sprunecriterion = 'mse';
else
    if (~isfield(sett, 'splitcriterion'))
        sett.splitcriterion = 'mse';
    end
    if (~isfield(sett, 'minleaf'))
        sett.minleaf = 1;
    end
    if (~isfield(sett, 'minparent'))
        sett.minparent = 10;
    end
    if (~isfield(sett, 'prune'))
        sett.prune = 'on';
    end
    if (~isfield(sett, 'prunecriterion'))
        sett.sprunecriterion = 'mse';
    end
end

if ((nargin < 5) || isempty(mdl))
    mdl_defined = false;
else
    mdl_defined = true;
end

%% Convert labels (if stored in cell array) into numbers
convert = false;

if (iscell(l_train))
    convert = true;
    unq_val = unique(l_train);
    l_train = conv_labels2mat(l_train);
end

%% Train the classifier
if (~mdl_defined)
    if verLessThan('matlab', '8.4')
        obj = RegressionTree.fit(d_train, l_train,                  ...
            'SplitCriterion', sett.splitcriterion,                  ...
            'MinLeaf', sett.minleaf,                                ...
            'MinParent', sett.minparent,                            ...
            'Prune', sett.prune,                                    ...
            'PruneCriterion', sett.prunecriterion);
    else
        obj = fitrtree(d_train, l_train,                            ...
            'SplitCriterion', sett.splitcriterion,                  ...
            'MinLeaf', sett.minleaf,                                ...
            'MinParent', sett.minparent,                            ...
            'Prune', sett.prune,                                    ...
            'PruneCriterion', sett.prunecriterion);
    end
else
    obj = mdl;
end
   
%% Evaluate the classifier
c_train = predict(obj, d_train);
c_test  = predict(obj, d_test);

%% Convert labels to proper representation (if needed)
if (convert)
    c_train = conv_mat2labels(c_train, unq_val);
    c_test  = conv_mat2labels(c_test, unq_val);
end