function score = check_classification_score(label)

% score = check_classifiation_score(label)
% 
% This function check chosen classification score string and returns 
% the method of classification scoring. 
%
% label     - scoring algorithm string (one of available possibilities)
% score     - output string containing the method of classif. scoring
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

%% Lower the score string (for consistency)
label = lower(label);

%% Set the classification specification (string)
if (strcmp(label, 'acc')                                            || ...
    strcmp(label, 'accuracy')                                       || ...
    strcmp(label, 'classification acc')                             || ...
    strcmp(label, 'classification accuracy'))

    % Classification accuracy
    score = 'acc';
    
elseif (strcmp(label, 'tss')                                        || ...
        strcmp(label, 'tradeoff sen spe')                           || ...
        strcmp(label, 'tradeoff between sen spe')                   || ...
        strcmp(label, 'tradeoff between sen and spe'))
    
    % Tradeoff between sensitivity and specificity
    score = 'tss';
    
elseif (strcmp(label, 'f1')                                         || ...
        strcmp(label, 'f1 label')                                   || ...
        strcmp(label, 'classification f1')                          || ...
        strcmp(label, 'classification f1 label'))
    
    % Classification F1 score
    score = 'f1';
    
elseif (strcmp(label, 'inf')                                        || ...
        strcmp(label, 'informedness')                               || ...
        strcmp(label, 'classification inf')                         || ...
        strcmp(label, 'classification informedness'))
    
    % Classification informedness
    score = 'inf';
    
elseif (strcmp(label, 'mar')                                        || ...
        strcmp(label, 'markedness')                                 || ...
        strcmp(label, 'classification mar')                         || ...
        strcmp(label, 'classification markedness'))
    
    % Classification markedness
    score = 'mar';
    
elseif (strcmp(label, 'pre')                                        || ...
        strcmp(label, 'prevalence')                                 || ...
        strcmp(label, 'classification pre')                         || ...
        strcmp(label, 'classification prevalence'))
    
    % Classification prevalence
    score = 'pre';
    
elseif (strcmp(label, 'sen')                                        || ...
        strcmp(label, 'sensitivity')                                || ...
        strcmp(label, 'classification sen')                         || ...
        strcmp(label, 'classification sensitivity'))
    
    % Classification sensitivity
    score = 'sen';
    
elseif (strcmp(label, 'spe')                                        || ...
        strcmp(label, 'specificity')                                || ...
        strcmp(label, 'classification spe')                         || ...
        strcmp(label, 'classification specificity'))
    
    % Classification specificity
    score = 'spe';
    
elseif (strcmp(label, 'dor')                                        || ...
        strcmp(label, 'odds ratio')                                 || ...
        strcmp(label, 'diagnostic odds')                            || ...
        strcmp(label, 'diagnostic odds ratio')                      || ...
        strcmp(label, 'classification odds')                        || ...
        strcmp(label, 'classification odds ratio'))
    
    % Classification (diagnostic) odds ratio
    score = 'dor';
    
elseif (strcmp(label, 'fdr')                                        || ...
        strcmp(label, 'false dr')                                   || ...
        strcmp(label, 'false discovery r')                          || ...
        strcmp(label, 'false discovery rate'))
    
    % False discovery rate
    score = 'fdr';
    
elseif (strcmp(label, 'fnr')                                        || ...
        strcmp(label, 'false nr')                                   || ...
        strcmp(label, 'false negative r')                           || ...
        strcmp(label, 'false negative rate'))
    
    % False negative rate
    score = 'fnr';
    
elseif (strcmp(label, 'for')                                        || ...
        strcmp(label, 'false or')                                   || ...
        strcmp(label, 'false omission r')                           || ...
        strcmp(label, 'false omission rate'))
    
    % False omission rate
    score = 'for';
    
elseif (strcmp(label, 'fpr')                                        || ...
        strcmp(label, 'false pr')                                   || ...
        strcmp(label, 'false positive r')                           || ...
        strcmp(label, 'false positive rate'))
    
    % False positive rate
    score = 'fpr';
    
elseif (strcmp(label, 'mcc')                                        || ...
        strcmp(label, 'mathews cc')                                 || ...
        strcmp(label, 'mathews corr coef')                          || ...
        strcmp(label, 'mathews correlation coef'))
    
    % Mathew's correlation coefficient
    score = 'mcc';
    
elseif (strcmp(label, 'nlr')                                        || ...
        strcmp(label, 'negative lr')                                || ...
        strcmp(label, 'negative likelihood r')                      || ...
        strcmp(label, 'negative likelihood ratio'))
    
    % Negative likelihood ration
    score = 'nlr';
    
elseif (strcmp(label, 'npv')                                        || ...
        strcmp(label, 'negative pv')                                || ...
        strcmp(label, 'negative predictive v')                      || ...
        strcmp(label, 'negative predictive value'))
    
    % Negative predictive value
    score = 'npv';
    
elseif (strcmp(label, 'plr')                                        || ...
        strcmp(label, 'positive lr')                                || ...
        strcmp(label, 'positive likelihood r')                      || ...
        strcmp(label, 'positive likelihood ratio'))
    
    % Positive likelihood ratio
    score = 'plr';
    
elseif (strcmp(label, 'ppv')                                        || ...
        strcmp(label, 'positive pv')                                || ...
        strcmp(label, 'positive predictive v')                      || ...
        strcmp(label, 'positive predictive value'))
    
    % Positive predictive value
    score = 'ppv';
    
elseif (strcmp(label, 'yjs')                                        || ...
        strcmp(label, 'youdens index')                              || ...
        strcmp(label, 'youdens statistics')                         || ...
        strcmp(label, 'youdens j statistics'))
    
    % Youden's J statistics
    score = 'yjs';
else
    score = [];
end