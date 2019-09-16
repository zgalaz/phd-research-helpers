%% test boxPlot
clear all;
close all;
clc;

addpath(genpath(pwd));

%% Lead the data
load('correlationPlot_testData.mat');

data = feat_matrix(:, 1);
clin = labels(:, 1);

%% Plot the correlation graph
options.xlabel = 'mean Noise-to-Harmonic Ratio';
options.ylabel = 'UPDRS III rating scale';

correlationPlot(data, cell2mat(clin), options);