function [c_train, c_test, obj] = perf_svmr(d_train, l_train, ...
    d_test, sett, mdl)

% [c_train, c_test, obj] = perf_svmr(d_train, l_train, d_test, sett)
% 
% This function performs regression (training; and testing) of the Support 
% Vector Machines (SVM) Regression algorithm. It returns predicted classes
% (multi-classes => prediction of 'continuous' values) for both training 
% and testing set of data. The classifier is specified by setings (sett) 
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
%   <http://www.mathworks.com/help/stats/fitrsvm.html>
%
% sett.kernel = Kernel function
%      {'linear', 'gaussian', 'rbf', 'polynomial'} 
%      (default: 'linear')
%
% Note: other settings will be implemented if necessary
%       (see the link above)
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
    sett.kernel = 'linear';
else
    if (~isfield(sett, 'kernel'))
        sett.kernel = 'linear';
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
    obj = fitrsvm(d_train, l_train, 'KernelFunction', sett.kernel);
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