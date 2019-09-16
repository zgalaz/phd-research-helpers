function [FDR] = false_discovery_rate(TP, FP, perc)

% [FDR] = false_discovery_rate(TP, FP, perc)
% 
% This function calculates the false discovery rate (FDR)
% TP        - true  positive
% FP        - false positive
% perc      - return the false discovery rate in percents
%             [0 = OFF, 1 = ON], default: 0 => FDR [-]
% FDR       - false discovery rate
%
% For more information see:
% http://en.wikipedia.org/wiki/False_discovery_rate
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

%% Paths and variables / calculate (FDR)
if ((nargin < 3) || (isempty(perc)) || perc == 0)
    FDR = (FP/(TP + FP));
else
    FDR = (FP/(TP + FP))*100;
end