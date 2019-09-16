%% Initial cleanup
clear all;
close all;
clc;

%% Paths and variables
addpath(genpath(pwd));

%% Load the *.mat file
load('test_reg.mat');

%% Preprocess the data
if (iscell(feat_matrix))
    feat_matrix = cell2mat(feat_matrix);
end
if (iscell(labels))
    labels = cell2mat(labels(:, 3:9));
end

feat_data = feat_matrix(:, 1:10);
clin_data = labels;

%% Perform the grid search
results = grid_search_cart(feat_data, clin_data);