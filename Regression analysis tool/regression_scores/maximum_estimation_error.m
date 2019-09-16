function [score] = maximum_estimation_error(actual, predic, max_num)

% [score] = maximum_estimation_error(actual, predic)
% 
% This function calculates the maximum estimation error (mean absolute 
% error / max. number in the vector)
%
% actual    - actual data
% predic    - predicted data
% max_num   - maximum number to divide MAE
%
% score     - maximum estimation error
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

%% Calculate the maximum estimation error
score = (mean_absolute_error(actual, predic)/max_num)*100;