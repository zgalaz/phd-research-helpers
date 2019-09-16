function gmm = GMMfit(data, labels, k)

% gmm = GMMfit(data, labels, k)
% 
% This function train Gaussian mixture model for each class in the data.
% 
% data          - input data matrix; columns are related to the features
%                 rows to the observations
% labels        - column vector with numeric labels of the observations
%                 (number of rows must eaqual to number of rows in data)
% k             - number of mixtures (default: 1)
% 
% gmm(i).class  - class ID
% gmm(i).obj    - object with the GMM model

%% Paths and variables
if ((nargin < 3) || isempty(k))
    k = 1;
end

%% Train the GMM for each class
classes = unique(labels);
gmm     = struct([]);

for i = 1:length(classes)
    gmm(i).class = classes(i);
    gmm(i).obj   = gmdistribution.fit(data(labels == classes(i), :), k);
end