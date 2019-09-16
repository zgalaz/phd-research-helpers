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
fun_sett.algorithm = 'bayes';
fun_sett.metrics   = {'acc', 'sen', 'spe'};
fun_sett.savetable = true;
fun_sett.tablename = 'learning_curves-bayes.xlsx';
fun_sett.plot      = true;

ml_sett.distrb = 'normal';
ml_sett.prior  = 'empirical';

%% Plot the learning curves
learning_curve_bayes( ...
    X(1:40, :),       ...
    X(41:end, :),     ...
    y(1:40),          ...
    y(41:end),        ...
    fun_sett,         ...
    ml_sett           ...
);