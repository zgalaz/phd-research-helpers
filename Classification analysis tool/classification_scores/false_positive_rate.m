function [FPR] = false_positive_rate(TN, FP, perc)

% [FPR] = false_positive_rate(TP, FP, perc)
% 
% This function calculates the false positive rate (FPR)
% TN        - true  negative
% FP        - false positive
% perc      - return the false positive rate in percents
%             [0 = OFF, 1 = ON], default: 0 => FPR [-]
% FPR       - false positive rate
%
% For more information see:
% http://en.wikipedia.org/wiki/False_positive_rate
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

%% Paths and variables / calculate (FPR)
if ((nargin < 3) || (isempty(perc)) || perc == 0)
    FPR = (FP/(TN + FP));
else
    FPR = (FP/(TN + FP))*100;
end