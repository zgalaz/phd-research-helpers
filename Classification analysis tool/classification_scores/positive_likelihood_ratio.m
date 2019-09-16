function [PLR] = positive_likelihood_ratio(TP, FN, perc)

% [PLR] = positive_likelihood_ratio(TP, FN, perc)
% 
% This function calculates the positive likelihood ratio (PLR)
% TP        - true  positive
% FN        - false negative
% perc      - return the positive likelihood ratio in percents
%             [0 = OFF, 1 = ON], default: 0 => PLR [-]
% PLR       - positive likelihood ratio
%
% For more information see:
% http://en.wikipedia.org/wiki/ ...
%   Likelihood_ratios_in_diagnostic_testing#positive_likelihood_ratio
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

%% Paths and variables / calculate (PLR)
TP = TP + 0.5;
FN = FN + 0.5;

if ((nargin < 3) || (isempty(perc)) || perc == 0)
    PLR  = (TP/FN);
else
    PLR  = (TP/FN)*100;
end