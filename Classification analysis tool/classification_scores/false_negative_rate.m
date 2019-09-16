function [FNR] = false_negative_rate(TP, FN, perc)

% [FNR] = false_negative_rate(TP, FN, perc)
% 
% This function calculates the false negative rate (FNR)
% TP        - true  positive
% FN        - false negative
% perc      - return the false negative rate in percents
%             [0 = OFF, 1 = ON], default: 0 => FNR [-]
% FNR       - false negative rate
%
% For more information see:
% http://en.wikipedia.org/wiki/ ...
%   False_positives_and_false_negatives#false_negative_rate
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

%% Paths and variables / calculate (FNR)
if ((nargin < 3) || (isempty(perc)) || perc == 0)
    FNR  = (FN/(TP + FN));
else
    FNR  = (FN/(TP + FN))*100;
end