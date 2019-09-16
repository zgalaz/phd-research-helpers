function [score] = mean_squared_error(actual, predic)

% [score] = mean_squared_error(actual, predic)
% 
% This function calculates the mean squared error
% actual    - actual data
% predic    - predicted data
%
% score     - mean squared error
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

%% Calculate the mean squared error
score = mean(squared_error(actual, predic));