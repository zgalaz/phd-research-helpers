%% Initial cleanup
clear all;
close all;
clc;

%% Paths and variables
addpath(genpath(pwd));

%% Load the *.mat file
load('test_data.mat');

y = round(rem((y + mean(y) + min(min(y))), 5)) + 1;

%% Prepare the structures
fun_sett = struct();
ml_sett  = struct();

fun_sett.normalize = true;
fun_sett.algorithm = 'or';
fun_sett.metrics   = {'mae', 'mse', 'rmse'};
fun_sett.savetable = true;
fun_sett.tablename = 'learning_curves-or.xlsx';
fun_sett.plot      = true;

ml_sett.link = 'logit';

%% Plot the learning curves
learning_curve_or( ...
    X(1:40, :),     ...
    X(41:end, :),   ...
    y(1:40),        ...
    y(41:end),      ...
    fun_sett,       ...
    ml_sett         ...
);