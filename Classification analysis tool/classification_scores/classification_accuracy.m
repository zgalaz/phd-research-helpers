function [ACC] = classification_accuracy(TP, FP, FN, TN, perc)

% [ACC] = classification_accuracy(TP, FP, FN, TN, perc)
% 
% This function calculates the classification accuracy (ACC)
% TP        - true  positive
% FP        - false positive
% FN        - false negative
% TN        - true  negative
% perc      - return the classification accuracy in percents
%             [0 = OFF, 1 = ON], default: 0 => ACC [-]
% ACC       - classification accuracy
%
% For more information see:
% http://en.wikipedia.org/wiki/Accuracy_and_precision
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

%% Paths and variables / calculate (ACC)
if ((nargin < 5) || (isempty(perc)) || perc == 0)
    ACC = ((TP + TN)/(TP + FP + FN + TN));
else
    ACC = ((TP + TN)/(TP + FP + FN + TN))*100;
end

