function [PPV] = positive_predictive_value(TP, FP, perc)

% [PPV] = positive_predictive_value(TP, FP, perc)
% 
% This function calculates the positive predictive value (PPV)
% TP        - true  positive
% FP        - false positive
% perc      - return the positive predictive value in percents
%             [0 = OFF, 1 = ON], default: 0 => PPV [-]
% PPV       - positive predictive value
%
% For more information see:
% http://en.wikipedia.org/wiki/Positive_and_negative_predictive_values
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

%% Paths and variables / calculate (PPV)
if ((nargin < 3) || (isempty(perc)) || perc == 0)
    PPV = (TP/(TP + FP));
else
    PPV = (TP/(TP + FP))*100;
end