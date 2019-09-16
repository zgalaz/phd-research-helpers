%% test clinicalCorr_polarPlot
clear all;
close all;
clc;

addpath(genpath(pwd));

%% 1. scenario: PD (9 clin. scales, 8 speech features, 1 measurement)
load('polarPlot_clinicalcorr_testData.mat');

fun = @clinicalCorr_polarPlot;
out = fun(cell2mat(speechFeatures(2:end, 1:5)), clinicalScales(:, 1:5));