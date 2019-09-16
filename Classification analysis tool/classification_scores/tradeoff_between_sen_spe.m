function [TSS] = tradeoff_between_sen_spe(TP, FP, FN, TN)

% [TSS] = tradeoff_between_sen_spe(TP, FP, FN, TN, perc)
% 
% This function calculates the classification tradeoff between 
% sensitivity and specificity (TSS)
% TP        - true  positive
% FP        - false positive
% FN        - false negative
% TN        - true  negative
% TSS       - tradeoff between sensitivity and specificity
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

%% Paths and variables / calculate (TSS)
x = classification_sensitivity(TP, FN, 0);
y = classification_specificity(FP, TN, 0);

TSS = 2^sin((pi*x)/2)*sin((pi*y)/2);
