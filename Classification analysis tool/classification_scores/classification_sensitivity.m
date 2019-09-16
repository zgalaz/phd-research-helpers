function [SEN] = classification_sensitivity(TP, FN, perc)

% [SEN] = classification_sensitivity(TP, FN, perc)
% 
% This function calculates the classification sensitivity (SEN)
% TP        - true  positive
% FN        - false negative
% perc      - return the classification sensitivity in percents
%             [0 = OFF, 1 = ON], default: 0 => SEN [-]
% SEN       - classification sensitivity
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

%% Paths and variables / calculate (SEN)
if ((nargin < 3) || (isempty(perc)) || perc == 0)
    SEN  = (TP/(TP + FN));
else
    SEN  = (TP/(TP + FN))*100;
end