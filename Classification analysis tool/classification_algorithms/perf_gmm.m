function [c_train, c_test, obj] = perf_gmm(d_train, l_train, ...
    d_test, sett, mdl)

% [c_train, c_test, obj] = perf_gmm(d_train, l_train, d_test, sett, mdl)
% 
% The functions performs classification (training and testing) of Gaussian
% Mixture Models classifier. This function returns predicted  classes for 
% both training and testing set of data. The classifier is  specified by
% setings (sett) structure.
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
%             the classification (testing step)
%
% OUTPUT DATA:
% c_train   - the classified training class (predicted labels)
% c_test    - the classified testing class (predicted labels)
% obj       - trained model object
% 
% CLASSIFIER OPTIONAL DATA:
% sett.mixt   = number of Gaussian mixtures 
%               (default: 1)
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
    sett.mixt = 1;
else
    if (~isfield(sett, 'mixt'))
        sett.mixt = 1;
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
    obj = GMMfit(d_train, l_train, sett.mixt);
else
    obj = mdl;
end   

%% Evaluate the classifier
c_train = GMMclassify(obj, d_train);
c_test  = GMMclassify(obj, d_test);

%% Convert labels to proper representation (if needed)
if (convert)
    c_train = conv_mat2labels(c_train, unq_val);
    c_test  = conv_mat2labels(c_test, unq_val);
end