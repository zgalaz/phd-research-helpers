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
fun_sett.algorithm = 'lda';
fun_sett.metrics   = {'acc', 'sen', 'spe'};
fun_sett.savetable = true;
fun_sett.tablename = 'learning_curves-lda.xlsx';
fun_sett.plot      = true;

ml_sett.discrm = 'linear';
ml_sett.prior  = 'empirical';

%% Plot the learning curves
learning_curve_lda( ...
    X(1:40, :),     ...
    X(41:end, :),   ...
    y(1:40),        ...
    y(41:end),      ...
    fun_sett,       ...
    ml_sett         ...
);