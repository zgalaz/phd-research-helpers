function [score] = mean_absolute_error(actual, predic)

% [score] = mean_absolute_error(actual, predic)
% 
% This function calculates the mean absolute error
% actual    - actual data
% predic    - predicted data
%
% score     - mean absolute error
%
%
%
% --
% ing. Zoltán Galáž
% xgalaz00@stud.feec.vutbr.cz       
% 
% Department of Telecommunications
% Faculty of Electrical Engineering and Communication
% Brno University of Technology

%% Calculate the mean absolute error
score = mean(absolute_error(actual, predic));