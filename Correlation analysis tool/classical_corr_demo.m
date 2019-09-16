%% Initial cleanup
clear all;
close all;
clc;

%% Paths and variables
addpath(genpath(pwd));

%% Load the *.mat file
load('test_corr.mat');

if (iscell(feat_matrix))
    feat_matrix = cell2mat(feat_matrix);
end

feat_labels = feat_labels(1:10).';
feat_data   = feat_matrix(:, 1:10);
clin_data   = clin_matrix;

%% Set the optional arguments
sett.type_corr = 'Spearman';

%% Perform correlation analysis
[rho, p_val] = correlation_analysis(feat_data, ...
    feat_labels, ...
    clin_matrix, ...
    sett);