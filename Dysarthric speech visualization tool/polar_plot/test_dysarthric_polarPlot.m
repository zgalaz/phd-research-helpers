%% test dysarthric_polarPlot
clear all;
close all;
clc;

addpath(genpath(pwd));

%% 1. scenario: PD/HC data (default settings, no feat. labels)
%     - 2  observations (1. row: PD; 2. row: HC)
%     - 16 prosodic features 
%           - 4 (Ph) features describing phonation
%           - 4 (Ar) features describing articulation
%           - 4 (Pr) features describing prosody
%           - 4 (Sr) features describing speech fluency
%     - 4  measurements in time
%
% -----------------------------------------------------------------------
% UNCOMMENT
dysarthria = struct();
dysarthria.Ph_features = rand(2, 6, 4);
dysarthria.Ar_features = rand(2, 5, 4);
dysarthria.Pr_features = rand(2, 9, 4);
dysarthria.Sr_features = rand(2, 4, 4);

out = dysarthric_polarPlot(dysarthria);

%% 2. scenario: multiple speakers data (default settings, no feat. labels)
%     - 2  observations (1. row: PD; 2. row: HC)
%     - 16 prosodic features
%           - 4 (Ph) features describing phonation
%           - 4 (Ar) features describing articulation
%           - 4 (Pr) features describing prosody
%           - 4 (Sr) features describing speech fluency
%     - 4  measurements in time
%
% -----------------------------------------------------------------------
% UNCOMMENT
% dysarthria = struct();
% dysarthria.Ph_features = rand(6, 4, 4); 
% dysarthria.Ar_features = rand(6, 4, 4);
% dysarthria.Pr_features = rand(6, 4, 4);
% dysarthria.Sr_features = rand(6, 4, 4);
% 
% out = dysarthric_polarPlot(dysarthria);