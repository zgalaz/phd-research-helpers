function [score] = root_mean_squared_log_error(actual, predic)

% [score] = root_mean_squared_log_error(actual, predic)
% 
% This function calculates the mean squared log error
% actual    - actual data
% predic    - predicted data
%
% score     - root mean squared log error
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

%% Calculate the root mean squared log error
score = sqrt(mean_squared_log_error(actual, predic));