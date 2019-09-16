function [score] = absolute_error(actual, predic)

% [score] = absolute_error(actual, predic)
% 
% This function calculates the absolute error
% actual    - actual data
% predic    - predicted data
%
% score     - absolute error
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

%% Calculate the absolute error
score = abs(actual(:) - predic(:));
