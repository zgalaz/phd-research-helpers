function [NLR] = negative_likelihood_ratio(FP, TN, perc)

% [NLR] = negative_likelihood_ratio(FP, TN, perc)
% 
% This function calculates the negative likelihood ratio (NLR)
% FP        - false  positive
% TN        - true negative
% perc      - return the negative likelihood ratio in percents
%             [0 = OFF, 1 = ON], default: 0 => NLR [-]
% NLR       - negative likelihood ratio
%
% For more information see:
% http://en.wikipedia.org/wiki/ ...
%   Likelihood_ratios_in_diagnostic_testing#negative_likelihood_ratio
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

%% Paths and variables / calculate (NLR)
FP = FP + 0.5;
TN = TN + 0.5;
if ((nargin < 3) || (isempty(perc)) || perc == 0)
    NLR  = (FP/TN);
else
    NLR  = (FP/TN)*100;
end