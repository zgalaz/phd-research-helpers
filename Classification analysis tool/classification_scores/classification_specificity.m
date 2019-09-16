function [SPE] = classification_specificity(FP, TN, perc)

% [SPE] = classification_specificity(FP, TN, perc)
% 
% This function calculates the classification specificity (SPE)
% FP        - false positive
% TN        - true negative
% perc      - return the classification specificity in percents
%             [0 = OFF, 1 = ON], default: 0 => SPE [-]
% SPE       - classification specificity
%
% For more information see:
% http://en.wikipedia.org/wiki/Sensitivity_and_specificity
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

%% Paths and variables / calculate (SPE)
if ((nargin < 3) || (isempty(perc)) || perc == 0)
    SPE  = (TN/(FP + TN));
else
    SPE  = (TN/(FP + TN))*100;
end

