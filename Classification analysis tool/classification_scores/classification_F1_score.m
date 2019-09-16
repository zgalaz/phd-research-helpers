function [F1] = classification_F1_score(TP, FP, FN, perc)

% [F1] = classification_F1_score(TP, FP, FN, perc)
% 
% This function calculates the F1 score (F1)
% TP        - true  positive
% FP        - false positive
% FN        - false negative
% perc      - return the classification F1 score in percents
%             [0 = OFF, 1 = ON], default: 0 => F1 [-]
% F1        - classification F1 score
%
% For more information see:
% http://en.wikipedia.org/wiki/F1_score
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

%% Paths and variables / calculate (F1)
if ((nargin < 4) || (isempty(perc)) || perc == 0)
    F1 = (2*TP/(2*TP + FP + FN));
else
    F1 = (2*TP/(2*TP + FP + FN))*100;
end