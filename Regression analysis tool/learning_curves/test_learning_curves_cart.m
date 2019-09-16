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
fun_sett.algorithm = 'cart';
fun_sett.metrics   = {'mae', 'mse', 'rmse'};
fun_sett.savetable = true;
fun_sett.tablename = 'learning_curves-tree.xlsx';
fun_sett.plot      = true;

ml_sett.splitcriterion  = 'mse';
ml_sett.minleaf         = 1;
ml_sett.minparent       = 10;
ml_sett.prune           = 'on';
ml_sett.sprunecriterion = 'error';

%% Plot the learning curves
learning_curve_cart( ...
    X(1:40, :),      ...
    X(41:end, :),    ...
    y(1:40),         ...
    y(41:end),       ...
    fun_sett,        ...
    ml_sett          ...
);