function [FOR] = false_omission_rate(TN, FN, perc)

% [FOR] = false_omission_rate(TN, FN, perc)
% 
% This function calculates the false omission rate (FOR)
% TN        - true  negative
% FN        - false negative
% perc      - return the positive false omission rate in percents
%             [0 = OFF, 1 = ON], default: 0 => FOR [-]
% FOM       - false omission rate
%
% For more information see:
% http://en.wikipedia.org/wiki/ ...
%   Positive_and_negative_predictive_values#false_omission_rate
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

%% Paths and variables / calculate (FOR)
if ((nargin < 3) || (isempty(perc)) || perc == 0)
    FOR = (FN/(TN + FN));
else
    FOR = (FN/(TN + FN))*100;
end