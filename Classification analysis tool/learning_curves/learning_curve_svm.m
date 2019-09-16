function [e_train, e_valid] = learning_curve_svm(X, X_val, y, y_val, ...
    fun_sett, ml_sett)

% [e_train, e_valid] = learning_curve_svm(X, X_val, y, y_val, ...
%    fun_sett, ml_sett)
% 
% This function trains/validates the support vector machine classification
% algorithm using predefined parameters in order to plot learning curves 
% (#training examples -> metrics). This enables us to detect if we face 
% the high bias or the high variance ML problem.
%
% INPUT:
% X         - input training data (rows: observations; cols: features)
% X_val     - input validation data (rows: observations; cols: features)
% y         - input training labels/dependent variable/s
% y_val     - input validation labels/dependent variable/s
% fun_sett  - input structure with the setting for selected algorithm 
% ml_sett   - input structure with the settings for learning process
%
% OUTPUT:
% e_train   - output training error series
% e_valid   - output validation error series
% 
% OPTIONAL DATA:
% fun_sett.normalize - data normalization (mean = 0, std = 1)
% fun_sett.algorithm - classification algorithm (predefined: svm)
% fun_sett.metrics   - metrics ecaluating the model ['acc', 'sen', etc.]
% fun_sett.savetable - save the result table to *.xlsx file
% fun_sett.tablename - name of the *.xlsx result table
% fun_sett.plot      - plot the results
% 
% CLASSIFIER OPTIONAL DATA:
% ml_sett.kernel = kernel used in SVM algorithm
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
if ((nargin < 3) || (isempty(fun_sett)))
    fun_sett.normalize = true;
    fun_sett.algorithm = 'svm';
    fun_sett.metrics   = {'acc', 'sen', 'spe'};
    fun_sett.savetable = true;
    fun_sett.tablename = 'learning_curves-svm.xlsx';
    fun_sett.plot      = true;
else
    if (~isfield(fun_sett, 'normalize'))
        fun_sett.normalize = true;
    end
    if (~isfield(fun_sett, 'algorithm'))
        fun_sett.algorithm = 'svm';
    end
    if (~isfield(fun_sett, 'metrics'))
        fun_sett.metrics = {'acc', 'sen', 'spe'};
    end
    if (~isfield(fun_sett, 'savetable'))
        fun_sett.savetable = true;
    end
    if (~isfield(fun_sett, 'tablename'))
        fun_sett.tablename = 'learning_curves-svm.xlsx';
    end
    if (~isfield(fun_sett, 'plot'))
        fun_sett.plot = true;
    end
end

if ((nargin < 4) || (isempty(ml_sett)))
    ml_sett.kernel = 'rbf';
else
    if (~isfield(ml_sett, 'kernel'))
        ml_sett.kernel = 'rbf';
    end
end

%% Preprocess the data
if (fun_sett.normalize)
    X     = zscore(X);
    X_val = zscore(X_val);
end

%% Randomly shuffle the data
perm_train = randperm(size(X, 1));
X = X(perm_train(:), :);
y = y(perm_train(:));

perm_valid = randperm(size(X_val, 1));
X_val = X_val(perm_valid(:), :);
y_val = y_val(perm_valid(:));

%% Divide the data (training/validation)
X_train = cell(size(X, 1), 1);
X_valid = cell(size(X, 1), 1);
y_train = cell(size(X, 1), 1);
y_valid = cell(size(X, 1), 1);

% Iterate over all train/validation splits
for i = 1:size(X, 1)
    
    % Select the data for each iteration
    X_train{i} = X(1:i, :);
    X_valid{i} = X_val;
    y_train{i} = y(1:i);
    y_valid{i} = y_val;
end

%% Compute the learning curve
metrics_train = nan(size(X, 1), length(fun_sett.metrics));
metrics_valid = nan(size(X, 1), length(fun_sett.metrics));

% Iterate over all train/validation splits
for i = 1:size(X, 1)
    
    % Define the classification settings
    cls.alg    = fun_sett.algorithm;
    cls.kernel = ml_sett.kernel;
                    
    % Train/validate the algorithm
    [lab_train, lab_valid, ~] = perf_classification( ...
        X_train{i}, ...
        y_train{i}, ...
        X_valid{i}, ...
        cls);
    
    % Compute the training confusion matrix
    true_D_train = (y_train{i} == 1);
    detc_D_train = (lab_train  == 1);
    true_H_train = (y_train{i} == 0);
    detc_H_train = (lab_train  == 0);

    TP_train = sum(true_D_train & detc_D_train);
    TN_train = sum(true_H_train & detc_H_train);
    FP_train = sum(true_H_train & detc_D_train);
    FN_train = sum(true_D_train & detc_H_train);
    
    % Compute the validation confusion matrix
    true_D_valid = (y_valid{i} == 1);
    detc_D_valid = (lab_valid  == 1);
    true_H_valid = (y_valid{i} == 0);
    detc_H_valid = (lab_valid  == 0);

    TP_valid = sum(true_D_valid & detc_D_valid);
    TN_valid = sum(true_H_valid & detc_H_valid);
    FP_valid = sum(true_H_valid & detc_D_valid);
    FN_valid = sum(true_D_valid & detc_H_valid);
    
    % Compute the classification metrics
    for m = 1:length(fun_sett.metrics)
        metrics_train(i, m) = calc_classification_score( ...
            TP_train, ...
            FP_train, ...
            FN_train, ...
            TN_train, ...
            fun_sett.metrics{m});
        metrics_valid(i, m) = calc_classification_score( ...
            TP_valid, ...
            FP_valid, ...
            FN_valid, ...
            TN_valid, ...
            fun_sett.metrics{m});
    end
end

% Set the output variables
e_train = metrics_train;
e_valid = metrics_valid;

%% Store the results to *.xlsx table
if (fun_sett.savetable)
    
    % Create the results table
    results_table = cell(size(X, 1) + 1, length(fun_sett.metrics)*2 + 2);
    
    % Prepare the training metrics header
    train_metrics = cell(1, length(fun_sett.metrics));
    for i = 1:length(fun_sett.metrics)
        train_metrics{i} = ['train (' fun_sett.metrics{i} ')'];
    end
    
    % Prepare the validation metrics header
    valid_metrics = cell(1, length(fun_sett.metrics));
    for i = 1:length(fun_sett.metrics)
        valid_metrics{i} = ['train (' fun_sett.metrics{i} ')'];
    end
    
    % Set the table header
    results_table(1, :) = [            ...
        {'training set'},              ...
        {'validation set'},            ...
        [train_metrics, valid_metrics] ...
    ];
    
    % Set the table values
    for i = 2:size(X, 1)      
        results_table(i, :) = [                ...
            num2cell(size(X_train{i - 1}, 1)), ...
            num2cell(size(X_valid{i - 1}, 1)), ...
            num2cell(e_train(i - 1, :)),       ...
            num2cell(e_valid(i - 1, :))        ...
        ];
    end
    
    % Save the table
    xlswrite(fun_sett.tablename, results_table);
end

%% Plot the results
if (fun_sett.plot)
    
    % Plot the training/validation metrics
    for i = 1:length(fun_sett.metrics)
        figure;
        
        hold on;
        grid on;
        plot(e_train(:, i), 'b');
        plot(e_valid(:, i), 'r');
        hold off;
        
        title(['SVM: learning curves (' fun_sett.metrics{i} ')']);
        xlabel('Number of training examples');
        ylabel(fun_sett.metrics{i});
        legend(['training ' fun_sett.metrics{i}], ...
               ['validation ' fun_sett.metrics{i}]);
    end
end