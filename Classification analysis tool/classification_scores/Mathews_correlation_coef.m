function [MCC] = Mathews_correlation_coef(TP, FP, FN, TN)

% [MCC] = Mathews_correlation_coef(TP, FP, FN, TN)
% 
% This function calculates the Mathews correlation coefficient (MCC)
% TP        - true  positive
% FP        - false positive
% FN        - false negative
% TN        - true  negative
% perc      - return the Mathews correlation coefficient in percents
%             [0 = OFF, 1 = ON], default: 0 => MCC [-]
% MCC       - Mathews correlation coefficient
%
% For more information see:
% http://en.wikipedia.org/wiki/Matthews_correlation_coefficient
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

%% Calculate the Matthew's corr. coeff.
nominator   = (TP*TN - FP*FN);
denominator = sqrt((TP + FP)*(TP + FN)*(TN + FP)*(TN + FN));

if (denominator == 0)
    denominator = 1;
end

MCC = nominator/denominator;
