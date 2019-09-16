function norm = normalize(data, method)

% norm = normalize(data, method)
% 
% This function normalizes the input data vector using several methods.
% Methods:
%   1) Standard Z-score normalization
%   2) MinMax normalization
%   3) SoftMax normalization
%   4) Sigmoid normalization
%
% data          - input data vector of values to normalize
% norm          - output data vector after normalization
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

%% Paths and variables
if ((nargin < 2) || isempty(method))
    method = 'Z-score';
end

%% Prepare the output data
norm = zeros(size(data, 1), size(data, 2));

%% Perform the normalization
if (strcmp(method, 'Z-score'))
    for dat = 1:size(data, 2)
        tmp = data(:, dat);
        norm(:, dat) = (tmp - mean(tmp))/std(tmp);
    end
elseif (strcmp(method, 'MinMax'))
    for dat = 1:size(data, 2)
        tmp = data(:, dat);
        norm(:, dat) = (tmp - min(tmp))/(max(tmp) - min(tmp));
    end
elseif (strcmp(method, 'SoftMax'))
    for dat = 1:size(data, 2)
        tmp = data(:, dat);
        den = 1 + exp(-((tmp - mean(tmp))/std(tmp)));
        norm(:, dat) = ones(size(data, 1), 1)./den;
    end
elseif (strcmp(method, 'Sigmoid'))
    for dat = 1:size(data, 2)
        tmp = data(:, dat);
        nom = 1 - exp(-((tmp - mean(tmp))/std(tmp)));
        den = 1 + exp(-((tmp - mean(tmp))/std(tmp)));
        norm(:, dat) = nom./den;
    end
else
    error('unsupported normalization method selected');
end
