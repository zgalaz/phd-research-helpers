%% Initial cleanup
clear all;
close all;
clc;

%% Paths and variables
addpath(genpath(pwd));

%% Load the *.mat file
load('test_cls.mat');

if (iscell(feat_matrix))
    feat_matrix = cell2mat(feat_matrix);
end
if (iscell(labels))
    temp = zeros(length(labels), 1);
    
    for i = 1:length(labels)
        if (strcmpi(labels{i}, 'pd'))
            temp(i, 1) = 1;
        else
            temp(i, 1) = 0;
        end
    end
    
    labels = temp;
    clear temp;
end

feat_data = feat_matrix(:, 1:10);
clin_data = labels;

%% Randomly shuffle the data
rand_ord  = randperm(length(clin_data));
feat_data = feat_data(rand_ord(:), :);
clin_data = clin_data(rand_ord(:));

%% Set the optional arguments
sett.cls.algorithm = 'tree';
sett.cls.metrics   = {'acc'};
sett.cv.validate   = true;
sett.cv.type_cv    = 'LeaveOut';
sett.cv.num_cv     = 1;
sett.cv.k_fold     = size(clin_data, 1);

%% Cross-validate the classification model
for label_idx  = 1:size(clin_data, 2)
    method.alg = sett.cls.algorithm;
    res_table  = cell(size(feat_data, 1) + 1, ...
        length(sett.cls.metrics)*2);
    
    for feature_idx = 1:size(feat_data, 2)
        disp(['Analysis : ',            ...
            ' (', num2str(feature_idx), ...
            '/', num2str(size(feat_data, 2)), ')']);
        
        % clin: training/testing labels
        % data: training/testing data
        clin = clin_data(:, label_idx);
        data = feat_data(:, feature_idx);

        % Perform the normalization (mean = 0, std = 1)
        norm = true;

        if (norm)
            for dat = 1:size(data, 2)   
                data_val = data(:, dat);
                mean_val = mean(data_val);
                std_val  = std(data_val);

                data(:, dat) = (data_val - mean_val)/std_val;
            end
        end

        % Cross-validate the selected classifier
        if (sett.cv.validate)
            
            % type_cv: type of cross-validation (k-fold/leave-on-out)
            % num_cv:  number of cross-validation runs
            % k_fold:  k-fold (if leave-one-out, k = size(data, 1))
            % p_labs:  predicted labels
            % p_labs:  true labels
            type_cv = sett.cv.type_cv;
            num_cv  = sett.cv.num_cv;
            k_fold  = sett.cv.k_fold;
            p_labs  = cell(num_cv, k_fold);
            t_labs  = cell(num_cv, k_fold);
        
            % Iterate over cross-validation runs
            for cv_idx = 1:num_cv
                predictions = cell(1, k_fold);
                presented   = cell(1, k_fold);

                % Get the cross-validation partitions
                if (strcmp(type_cv, 'KFold'))
                    cv_paritions = ...
                        cvpartition(numel(clin), 'KFold', k_fold);
                else
                    cv_paritions = ...
                        cvpartition(numel(clin), 'LeaveOut');
                end

                % Train and evaluate the classifier
                for fold_idx  = 1:k_fold
                    train_idx = cv_paritions.training(fold_idx);
                    test_idx  = cv_paritions.test(fold_idx);

                    train_table = data(train_idx);
                    train_label = clin(train_idx);
                    test_table  = data(test_idx);
                    test_label  = clin(test_idx);

                    [~, pred_nested, cls_model] = ...
                        perf_classification( ...
                            train_table, ...
                            train_label, ...
                            test_table,  ...
                            method);

                    predictions{fold_idx} = pred_nested(:);
                    presented{fold_idx}   = test_label(:);
                end
                
                p_labs(cv_idx, :) = predictions;
                t_labs(cv_idx, :) = presented;
            end
            
            % Compute the metrics
            for metric_idx = 1:length(sett.cls.metrics)
                metric = zeros(num_cv, k_fold);
                
                for cv_idx = 1:num_cv
                    for fold_idx  = 1:k_fold
                        true_labs = t_labs{cv_idx, fold_idx};
                        pred_labs = p_labs{cv_idx, fold_idx};
                        
                        true_D = (t_labs{cv_idx, fold_idx} == 1);
                        detc_D = (p_labs{cv_idx, fold_idx} == 1);
                        true_H = (t_labs{cv_idx, fold_idx} == 0);
                        detc_H = (p_labs{cv_idx, fold_idx} == 0);

                        TP = sum(true_D & detc_D);
                        TN = sum(true_H & detc_H);
                        FP = sum(true_H & detc_D);
                        FN = sum(true_D & detc_H);
                    
                        metric(cv_idx, fold_idx) = ...
                            calc_classifiation_score( ...
                                TP, ...
                                FP, ...
                                FN, ...
                                TN, ...
                                sett.cls.metrics{metric_idx});
                    end
                end
                
                % Set the result values
                res_table{feature_idx + 1, ...
                    (metric_idx - 1)*2 + 1} = mean(metric(:));
                res_table{feature_idx + 1, ...
                    (metric_idx - 1)*2 + 2} = std(metric(:));
                disp([' ', sett.cls.metrics{metric_idx}, ... 
                    ' = ', num2str(mean(metric(:))), ...
                    '+-', num2str(std(metric(:)))]);
            end
        else
            error('TODO: no cross-validation');
        end
    end
    
    xlswrite(['analysis-' num2str(label_idx) '.xlsx'], res_table);
end