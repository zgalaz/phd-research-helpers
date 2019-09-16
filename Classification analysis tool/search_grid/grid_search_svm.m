function res_struct = grid_search_svm(data, clin, sett, search)

% res_struct = grid_search_svm(data, clin, sett, search)
% 
% This function performs a grid search on support vector machines (svm) 
% classifier using predefined parameters. The results of the function 
% are stored into the output structure and can also be stored in *.xlsx 
% file.
%
% INPUT DATA:
% data           - input matrix [columns: features; rows: observations]
% clin           - input vector of labels (length(clin) == size(data, 1))
% sett           - structure with the setting for selected algorithm
% search         - structure with the settings for grid search
%
% OUTPUT DATA:
% res_struct     - output structure with the results
%
% OPTIONAL DATA:
% sett.normalize - data normalization (mean = 0, std = 1)
% sett.algorithm - classification algorithm (predefines: svm)
% sett.metrics   - metrics ecaluating the model ['acc', 'sen', etc.]
% sett.type_cv   - type of cross validation {'KFold' 'LeaveOut'}
% sett.num_cv    - number of cross validation runs
% sett.k_fold    - k-fold validation
% sett.savetable - save the result table to *.xlsx file
% sett.tablename - name of the *.xlsx result table
% sett.disp      - display the results to console
%
% CLASSIFIER OPTIONAL DATA:
% search.kernel = kernel used in SVM algorithm
%               {'linear', 'polynomial', 'rbf'}
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
if ((nargin < 3) || (isempty(sett)))
    sett.normalize = true;
    sett.algorithm = 'svm';
    sett.metrics   = {'mcc', 'acc', 'sen', 'spe'};
    sett.type_cv   = 'KFold';
    sett.num_cv    = 10;
    sett.k_fold    = 10;
    sett.savetable = true;
    sett.tablename = 'grid_search-svm.xlsx';
    sett.disp      = true;
else
    if (~isfield(sett, 'normalize'))
        sett.normalize = true;
    end
    if (~isfield(sett, 'algorithm'))
        sett.algorithm = 'svm';
    end
    if (~isfield(sett, 'metrics'))
        sett.metrics = {'mcc', 'acc', 'sen', 'spe'};
    end
    if (~isfield(sett, 'type_cv'))
        sett.type_cv = 'KFold';
    end
    if (~isfield(sett, 'num_cv'))
        sett.num_cv = 10;
    end
    if (~isfield(sett, 'k_fold'))
        sett.k_fold = 10;
    end
    if (~isfield(sett, 'savetable'))
        sett.savetable = true;
    end
    if (~isfield(sett, 'tablename'))
        sett.tablename = 'grid_search-svm.xlsx';
    end
    if (~isfield(sett, 'disp'))
        sett.disp = true;
    end
end

if ((nargin < 4) || (isempty(search)))
    search.kernel = {'linear', 'polynomial', 'rbf'};
else
    if (~isfield(search, 'kernel'))
        search.kernel = {'linear', 'polynomial', 'rbf'};
    end
end

%% Preprocess the data
% Convert data/labels to numerical form
if (iscell(data))
    data = cell2mat(data);
end
if (iscell(clin))
    clin = cell2mat(clin);
end

% Clean the data/labels
data  = data';
nans  = any(isnan(data));
infs  = any(isinf(data));
comp  = any(imag(data) > 0);
zero  = sum(abs(data)) == 0;
leave = ~(nans | infs | comp | zero) == 1;

data  = data(:, leave);
data  = data';
clin  = clin(leave);

% Normalize the data
if (sett.normalize)
    for dat = 1:size(data, 2)   
        data_val = data(:, dat);
        mean_val = mean(data_val);
        std_val  = std(data_val);

        data(:, dat) = (data_val - mean_val)/std_val;
    end
end

% Randomly shuffle the data/labels
rord = randperm(length(clin));
data = data(rord(:), :);
clin = clin(rord(:));

%% Set the grid search parameters
sett.g_search = {'kernel', search.kernel};

%% Get the number of combinations in the grid
p1 = 1:length(sett.g_search{1, 2});

dimensions = cellfun(@numel, {p1});
[i1] = ind2sub(dimensions, 1:prod(dimensions));
combinations = size([p1(i1)]', 1);

%% Prepare the output structure/table
res_struct = struct();
res_table  = cell(combinations + 1, length(sett.metrics)*2);
descript   = cell(1, combinations);
res_iter   = 1;

%% Perform the grid search
for p1 = 1:length(sett.g_search{1, 2})
        
    % Set the temporary settings
    cls.alg    = sett.algorithm;
    cls.kernel = sett.g_search{1, 2}{p1};

    % Create the display string
    if (sett.disp)
        descript{res_iter} = [sett.g_search{1, 1}, '(', cls.kernel, ')'];
        
        display  = ['Classifier (svm), ', ...
            sett.g_search{1, 1}, '(', cls.kernel, ') -> ', ...
        ];
    end

    % type_cv: type of cross-validation (k-fold/leave-on-out)
    % num_cv:  number of cross-validation runs
    % k_fold:  k-fold (if leave-one-out, k = size(data, 1))
    % p_labs:  predicted labels
    % p_labs:  true labels
    type_cv = sett.type_cv;
    num_cv  = sett.num_cv;
    k_fold  = sett.k_fold;
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
                    cls);

            predictions{fold_idx} = pred_nested(:);
            presented{fold_idx}   = test_label(:);
        end

        p_labs(cv_idx, :) = predictions;
        t_labs(cv_idx, :) = presented;
    end

    % Compute the metrics
    for metric_idx = 1:length(sett.metrics)
        metric  = zeros(num_cv, k_fold);
        h_means = (metric_idx - 1)*2 + 1 + length(sett.g_search);
        h_stds  = (metric_idx - 1)*2 + 2 + length(sett.g_search);

        for cv_idx = 1:num_cv
            for fold_idx  = 1:k_fold
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
                        sett.metrics{metric_idx});
            end
        end

        % Post-process the metrics
        metric = metric(:);
        metric = metric(~isnan(metric));

        metric_mean = mean(metric);
        metric_std  = std(metric);

        % Update the results table
        res_table{1, 1} = sett.g_search{1, 1};
        res_table{1, h_means} = ...
            [sett.metrics{metric_idx} ' (mean)'];
        res_table{1, h_stds} = ...
            [sett.metrics{metric_idx} ' (std)'];
        res_table{res_iter + 1, 1} = cls.kernel;
        res_table{res_iter + 1, h_means} = metric_mean;
        res_table{res_iter + 1, h_stds}  = metric_std;

        % Update the results structure
        res_struct(res_iter).res_table = res_table;
        res_struct(res_iter).cls_model = cls_model;

        % Update the display string
        if (sett.disp)
            display = [display, ...
                ' ', sett.metrics{metric_idx}, ... 
                ' = ', num2str(mean(metric(:))), ...
                '+-', num2str(std(metric(:))), '; ', ...
            ];
        end
    end

    % Display the results
    if (sett.disp)
        disp(display);
    end

    % Update the iterator
    res_iter = +res_iter + 1;
end

%% Plot the errorbar graphs
for metric_idx = 1:length(sett.metrics)
    h_means = (metric_idx - 1)*2 + 1 + length(sett.g_search);
    h_stds  = (metric_idx - 1)*2 + 2 + length(sett.g_search);
    d_means = cell2mat(res_table(2:end, h_means));
    d_stds  = cell2mat(res_table(2:end, h_stds));

    figure;
    errorbar(d_means, d_stds);
    
    ticks = combinations;
    lims  = get(gca, 'XLim');
    set(gca(), 'XTick', linspace(lims(1), lims(2), ticks))
    
    set(gca(), 'XTickLabel', descript);
    rotateXLabels(gca(), 45);
    
    title(['Classifier (svm): ' ...
        sett.metrics{metric_idx} ' (mean +- std)']);
end

%% Save the results
if (sett.savetable)
    xlswrite(sett.tablename, res_table);
end