function [c_train, c_test, obj] = perf_or(d_train, l_train, ...
    d_test, sett, mdl)

% [c_train, c_test, obj] = perf_or(d_train, l_train, d_test, sett)
% 
% The functions performs the regression task (training; testing) of Ordinal
% Regression (OR) algorithm, and returns the predicted classes' probability
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
% c_train   - the classified training class probability
% c_test    - the classified testing class probability
% obj       - trained model object
%
% CLASSIFIER OPTIONAL DATA:
% for more information, see: 
%   <http://www.mathworks.com/help/stats/mnrfit.html>
%
% sett.model = Type of model to fit 
%              {'nominal' 'ordinal' 'hierarchical'}
%              (default: 'ordinal')
% sett.link  = Link function 
%              {'logit' 'probit' 'comploglog' 'loglog'}
%              (default: 'logit')
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

%% Paths and variables
if ((nargin < 4) || (isempty(sett)))
    sett.model = 'ordinal';
    sett.link  = 'logit';
else
    if (~isfield(sett, 'model'))
        sett.model = 'ordinal';
    end
    if (~isfield(sett, 'link'))
        
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
    [B, dev, stats] = mnrfit(d_train, l_train,                      ...
        'Model', sett.model,                                        ...
        'Link', sett.link);
    
    obj.stats = stats;
    obj.dev   = dev;
    obj.B     = B;
end

%% Evaluate the classifier
[pihat_train, ~, ~] = mnrval(B, d_train, stats, 'Model', sett.model);
[pihat_test, ~, ~]  = mnrval(B, d_test, stats, 'Model', sett.model);

[~, ind_train] = max(pihat_train');
[~, ind_test]  = max(pihat_test');

c_train = ind_train(:);
c_test  = ind_test(:);

%% Convert labels to proper representation (if needed)
if (convert)
    c_train = conv_mat2labels(c_train, unq_val);
    c_test  = conv_mat2labels(c_test, unq_val);
end