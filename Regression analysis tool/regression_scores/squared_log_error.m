function [score] = squared_log_error(actual, predic)

% [score] = squared_log_error(actual, predic)
% 
% This function calculates the squared log error
% actual    - actual data
% predic    - predicted data
%
% score     - squared log error
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

%% Calculate the squared log error
score = (log(1 + actual(:)) - log(1 + predic(:))).^2;