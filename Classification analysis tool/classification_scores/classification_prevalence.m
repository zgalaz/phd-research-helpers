function [PRE] = classification_prevalence(TP, FP, FN, TN, perc)

% [PRE] = classification_prevalence(TP, FP, FN, TN, perc)
% 
% This function calculates the classification prevalence (PRE)
% TP        - true  positive
% FP        - false positive
% FN        - false negative
% TN        - true  negative
% perc      - return the classification prevalence in percents
%             [0 = OFF, 1 = ON], default: 0 => PRE [-]
% ACC       - classification prevalence
%
% For more information see:
% http://en.wikipedia.org/wiki/Prevalence
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

%% Paths and variables / calculate (PRE)
if ((nargin < 5) || (isempty(perc)) || perc == 0)
    PRE = ((TP + FN)/(TP + FP + FN + TN));
else
    PRE = ((TP + FN)/(TP + FP + FN + TN))*100;
end

