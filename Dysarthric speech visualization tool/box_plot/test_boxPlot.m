%% test boxPlot
clear all;
close all;
clc;

addpath(genpath(pwd));

%% 1. scenario
%     - multiple (100) data values
%     - multiple (2) groups
%
% -----------------------------------------------------------------------
% UNCOMMENT 
% data = rand(100, 2);
% 
% boxPlot(data);

%% 2. scenario
%     - multiple (100) data values
%     - multiple (10) groups
%
% -----------------------------------------------------------------------
% UNCOMMENT 
% data = rand(100, 10);
% 
% boxPlot(data);

%% 3. scenario (web-scenario)
load('boxPlot_testData.mat');

idx  = 14:18;
data = feat_matrix(:, idx);

options.title  = 'Box plot (speech features) for PD patients';
options.xlabel = 'Speech features';
options.ylabel = 'Feature values';
options.labels = {'NVB',                ...
                  'DVB',                ...
                  'TST',                ...
                  'NST',                ...
                  'TPT'};

boxPlot(data, options);