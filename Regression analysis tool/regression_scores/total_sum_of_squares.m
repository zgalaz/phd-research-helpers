function [score] = total_sum_of_squares(actual)

% [score] = total_sum_of_squares(actual)
% 
% This function calculates the total sum of squares.
% actual    - actual data
%
% score     - total sum of squares
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

%% Calculate the total sum of squares
score = sum((actual(:) - mean(actual(:))).^2);