%% test general_polarPlot
clear all;
close all;
clc;

addpath(genpath(pwd));

%% 1. scenario: PD/HC data (default settings, no feat. labels)
%     - 2  observations (1. row: PD; 2. row: HC)
%     - 16 speech features 
%     - 4  measurements in time
%
% -----------------------------------------------------------------------
% UNCOMMENT 
data = rand(2, 16, 1);
out  = general_polarPlot(data);

%% 1. scenario: multiple speakers data (default settings, no feat. labels)
%     - 6  observations
%     - 16 speech features 
%     - 4  measurements in time
%
% -----------------------------------------------------------------------
% UNCOMMENT 
% data = rand(6, 16, 4);
% out  = general_polarPlot(data);