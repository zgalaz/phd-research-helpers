%% test boxPlot
clear all;
close all;
clc;

addpath(genpath(pwd));

%% Lead the data
load('kernel_densityPlot_testData.mat');

options.xlabel = 'pitch period entropy';
options.ylabel = 'probability density function';
options.legend = {'healthy controls', 'PD patients'};
options.title  = 'Kernel density estimation';

data    = cell(1, 2);
data{1} = feat_matrix(1:27, 3);
data{2} = feat_matrix(28:end, 3);

%% Plot the kernel probability density estimation graph
kernel_densityPlot(data, options);