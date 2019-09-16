function [DOR, SE, p] = diagnostic_odds_ratio(TP, FP, FN, TN, logarithm)

% [DOR, SE, p] = diagnostic_odds_ratio(TP, FP, FN, TN, logarithm)
% 
% This function calculates the diagnostic odds ratio (DOR)
% TP        - true  positive
% FP        - false positive
% FN        - false negative
% TN        - true  negative
% perc      - return the diagnostic odds ratio in percents
%             [0 = OFF, 1 = ON], default: 0 => DOR [-]
% DOR       - diagnostic odds ratio
% SE        - standard error
% p         - 95% confidence interval
%
% For more information see:
% http://en.wikipedia.org/wiki/Diagnostic_odds_ratio
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

%% Paths and variables
p  = zeros(1, 2);
SE = 0;

%% Avoid: Criticism (from Wiki)
%         "diagnostic odds ratio is undefined when the number
%          of false negatives or false positives is zero"
TP = TP + 0.5;
TN = TN + 0.5;
FP = FP + 0.5;
FN = FN + 0.5;

%% Calculate (DOR)
if ((nargin < 5) || (isempty(logarithm)) || logarithm == 0)
    DOR  = (TP/FN)/(FP/TN);
else
    DOR  = log((TP/FN)/(FP/TN));
    SE   = sqrt(1/TP + 1/FP + 1/FN + 1/TN);
    p(1) = DOR - 1.96*SE;
    p(2) = DOR + 1.96*SE;
end
