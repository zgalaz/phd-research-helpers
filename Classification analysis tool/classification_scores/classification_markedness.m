function [MAR] = classification_markedness(TP, FP, FN, TN, perc)

% [MAR] = classification_markedness(TP, FP, FN, TN, perc)
% 
% This function calculates the classification markedness (MAR)
% TP        - true  positive
% FP        - false positive
% FN        - false negative
% TN        - true  negative
% perc      - return the classification markedness in percents
%             [0 = OFF, 1 = ON], default: 0 => MAR [-]
% MAR       - classification classification markedness
%
% For more information see:
% http://en.wikipedia.org/wiki/Markedness
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
if ((nargin < 5) || (isempty(perc)) || perc == 0)
    PPV = positive_predictive_value(TP, FP, 0);
    NPV = negative_predictive_value(TN, FN, 0);

else
    PPV = positive_predictive_value(TP, FP, 1);
    NPV = negative_predictive_value(TN, FN, 1);
end

%% Calculate (MAR)
MAR = PPV + NPV - 1;