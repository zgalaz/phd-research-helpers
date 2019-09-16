function [score] = coefficient_of_determination(actual, predic)

% [score] = coefficient_of_determination(actual, predic)
% 
% This function calculates the  r squared.
% actual    - actual data
% predic    - predicted data
%
% score     - r squared
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

%% Calculate the coefficient_of_determination
num = residual_sum_of_squares(actual, predic);
den = total_sum_of_squares(actual);

score = 1 - num/den;