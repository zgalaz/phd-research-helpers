function [score] = ...
    relative_estimation_error(actual, predic, data_range, perc)

% [score] = relative_estimation_error(actual, predic, data_range, perc)
% 
% This function calculates the relative estimation error. It is computed
% as: mean squared error/range of data.
%
% actual        - actual data
% predic        - predicted data
% data_range    - range of data
% perc          - return the score in percents
%                 [0 = OFF, 1 = ON], default: 0 => score [-]
%
% score         - squared error
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

%% Paths and variables / calculate (EER)
if ((nargin < 4) || (isempty(perc)) || perc == 0)
    score = mean_absolute_error(actual, predic)/data_range;
else
    score = mean_absolute_error(actual, predic)/data_range*100;
end
