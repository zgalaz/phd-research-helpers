function [YJS] = Youdens_J_statistics(TP, FP, FN, TN)

% [YJS] = Youdens_J_statistics(TP, FP, FN, TN)
% 
% This function calculates the Youden's J statistics (YJS)
% TP        - true  positive
% FP        - false positive
% FN        - false negative
% TN        - true  negative
% perc      - return the Mathews correlation coefficient in percents
%             [0 = OFF, 1 = ON], default: 0 => MCC [-]
% MCC       - Youden's J statistics
%
% For more information see:
% https://en.wikipedia.org/wiki/Youden%27s_J_statistic
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

%% Paths and variables / calculate (MCC)
SEN = classification_sensitivity(TP, FN, 0);
SPE = classification_specificity(FP, TN, 0);
YJS = SEN + SPE - 1;
