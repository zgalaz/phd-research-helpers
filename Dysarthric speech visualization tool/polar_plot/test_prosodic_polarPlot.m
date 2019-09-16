%% test dysarthric_polarPlot
clear all;
close all;
clc;

addpath(genpath(pwd));

%% 1. scenario: PD/HC data (default settings, no feat. labels)
%     - 2  observations (1. row: PD; 2. row: HC)
%     - 12 prosodic features 
%           - 4 (F0) features describing monopitch
%           - 4 (A0) features describing monoloudness
%           - 4 (SR) features describing speech rate
%     - 4  measurements in time
%
% -----------------------------------------------------------------------
% UNCOMMENT
% prosody = struct();
% prosody.F0_features = rand(2, 4, 4); 
% prosody.A0_features = rand(2, 4, 4);
% prosody.SR_features = rand(2, 4, 4);
% 
% out = prosodic_polarPlot(prosody);

%% 2. scenario: multiple speakers data (default settings, no feat. labels)
%     - 3  observations
%     - 12 prosodic features 
%           - 4 (F0) features describing monopitch
%           - 4 (A0) features describing monoloudness
%           - 4 (SR) features describing speech rate
%     - 4  measurements in time
%
% -----------------------------------------------------------------------
% UNCOMMENT
prosody = struct();
prosody.F0_features = rand(3, 4, 4); 
prosody.A0_features = rand(3, 4, 4);
prosody.SR_features = rand(3, 4, 4);

out = prosodic_polarPlot(prosody);