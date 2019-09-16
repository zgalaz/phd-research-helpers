%% Initial cleanup
clear all;
close all;
clc;

%% Paths and variables
addpath(genpath(pwd));

%% Load the *.mat file
load('test_cls.mat');

%% Preprocess the data
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

%% Perform the grid search
results = grid_search_knn(feat_data, clin_data);