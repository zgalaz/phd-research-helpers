%% Initial cleanup
clear all;
close all;
clc;

%% Paths and variables
addpath(genpath(pwd));

%% Load the *.mat file
load('test_data.mat');

%% Prepare the structures
fun_sett = struct();
ml_sett  = struct();

fun_sett.normalize = true;
fun_sett.algorithm = 'svmr';
fun_sett.metrics   = {'mae', 'mse', 'rmse'};
fun_sett.savetable = true;
fun_sett.tablename = 'learning_curves-svmr.xlsx';
fun_sett.plot      = true;

ml_sett.kernel = 'gaussian';

%% Plot the learning curves
learning_curve_svmr( ...
    X(1:40, :),     ...
    X(41:end, :),   ...
    y(1:40),        ...
    y(41:end),      ...
    fun_sett,       ...
    ml_sett         ...
);