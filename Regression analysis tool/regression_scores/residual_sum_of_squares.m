function [score] = residual_sum_of_squares(actual, predic)

% [score] = residual_sum_of_squares(actual, predic)
% 
% This function calculates the residual sum of squares.
% actual    - actual data
% predic    - predicted data
%
% score     - residual sum of squares
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

%% Calculate the residual sum of squares
score = sum((actual(:) - predic(:)).^2);