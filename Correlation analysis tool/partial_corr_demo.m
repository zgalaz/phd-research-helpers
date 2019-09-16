%% Initial cleanup
clear all;
close all;
clc;

%% Paths and variables
addpath(genpath(pwd));

%% Load the *.mat file
load('test_data.mat');

%% Prepare the data
feature_matrix  = {feat_matrix_1, feat_matrix_2};
clinical_matrix = {clin_matrix_1, clin_matrix_2};
covariates      = {covs_matrix_1, covs_matrix_2};

%% Set the optional arguments
options.type_corr    = {'Pearson', 'Spearman'};
options.analyse_BFBC = true;
options.analyse_BFDC = true;
options.analyse_DFDC = true;
options.feat_labels  = feat_labels;
options.clin_labels  = clin_labels;
options.create_xlsx  = true;

%% Perform correlation analysis
out_results = partial_correlation( ...
    feature_matrix,                ...
    clinical_matrix,               ...
    covariates,                    ...
    options                        ...
);

%% Plot residuals/diagnostics (example)
plotResiduals(out_results.lin_models{1, 1}.feat_data_mdl(1).lin_mdl);
plotDiagnostics(out_results.lin_models{1, 1}.feat_data_mdl(1).lin_mdl);