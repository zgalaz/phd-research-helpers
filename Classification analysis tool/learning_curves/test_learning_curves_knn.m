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
fun_sett.algorithm = 'knn';
fun_sett.metrics   = {'acc', 'sen', 'spe'};
fun_sett.savetable = true;
fun_sett.tablename = 'learning_curves-knn.xlsx';
fun_sett.plot      = true;

ml_sett.k    = 2;
ml_sett.dist = 'euclidean';

%% Plot the learning curves
learning_curve_knn( ...
    X(1:40, :),     ...
    X(41:end, :),   ...
    y(1:40),        ...
    y(41:end),      ...
    fun_sett,       ...
    ml_sett         ...
);