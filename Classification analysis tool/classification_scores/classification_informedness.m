function [INF] = classification_informedness(TP, FP, FN, TN, perc)

% [INF] = classification_informedness(TP, FP, FN, TN, perc)
% 
% This function calculates the classification informedness (INF)
% TP        - true  positive
% FP        - false positive
% FN        - false negative
% TN        - true  negative
% perc      - return the classification informedness in percents
%             [0 = OFF, 1 = ON], default: 0 => INF [-]
% INF       - classification classification informedness
%
% For more information see:
% http://en.wiktionary.org/wiki/informedness
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
    SEN = classification_sensitivity(TP, FN, 0);
    SPE = classification_specificity(FP, TN, 0);

else
    SEN = classification_sensitivity(TP, FN, 1);
    SPE = classification_specificity(FP, TN, 1);
end

%% Calculate (INF)
INF = SEN + SPE - 1;