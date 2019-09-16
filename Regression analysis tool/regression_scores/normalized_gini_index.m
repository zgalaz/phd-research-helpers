function [score] = normalized_gini_index(actual, predic)

% [score] = normalized_gini_index(actual, predic)
%
% This function calculates the gini index
% actual    - actual data (n*1 matrix of actual values)
% predic    - predicted data (n*1 matrix of predicted values)
%
% score     - normalized gini index
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

%% Calculate the normalized gini index
score = gini_index(actual, predic)/gini_index(actual, actual);