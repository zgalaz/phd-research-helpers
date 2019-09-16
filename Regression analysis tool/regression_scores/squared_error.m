function [score] = squared_error(actual, predic)

% [score] = squared_error(actual, predic)
% 
% This function calculates the squared error
% actual    - actual data
% predic    - predicted data
%
% score     - squared error
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

%% Calculate the squared error
score = (actual(:) - predic(:)).^2;