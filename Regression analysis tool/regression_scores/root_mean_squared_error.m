function [score] = root_mean_squared_error(actual, predic)

% [score] = root_mean_squared_error(actual, predic)
% 
% This function calculates the root mean squared error
% actual    - actual data
% predic    - predicted data
%
% score     - root mean squared error
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

%% Calculate the root mean squared error
score = sqrt(mean_squared_error(actual, predic));