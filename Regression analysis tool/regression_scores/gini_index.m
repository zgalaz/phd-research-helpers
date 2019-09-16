function [score] = gini_index(actual, predic)

% [score] = gini_index(actual, predic)
%
% This function calculates the gini index
% actual    - actual data (n*1 matrix of actual values)
% predic    - predicted data (n*1 matrix of predicted values)
%
% score     - gini index
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

%% Sort predicted values (multiply by (-1))
[~, I] = sort(-predic);

%% Input variables
delta = 1/length(actual);
APPS  = 0;
ALPS  = 0;
score = 0;

%% Prepare variable to hold total number of losses
total_losses = sum(actual);

%% Calculate the gini index
for i = 1:length(actual)
    loc = I(i);
    
    ALPS = ALPS + actual(loc)/total_losses;
    APPS = APPS + delta;
    
    score = score + ALPS - APPS;
end

score = score/length(actual);