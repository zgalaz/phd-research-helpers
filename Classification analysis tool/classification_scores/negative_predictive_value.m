function [NPV] = negative_predictive_value(TN, FN, perc)

% [NPV] = negative_predictive_value(TN, FN, perc)
% 
% This function calculates the negative predictive value (NPV)
% FN        - false negative
% TN        - true  negative
% perc      - return the negative predictive value in percents
%             [0 = OFF, 1 = ON], default: 0 => NPV [-]
% NPV       - negative predictive value
%
% For more information see:
% http://en.wikipedia.org/wiki/Positive_and_negative_predictive_values
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

%% Paths and variables / calculate (NPV)
if ((nargin < 3) || (isempty(perc)) || perc == 0)
    NPV = (TN/(TN + FN));
else
    NPV = (TN/(TN + FN))*100;
end