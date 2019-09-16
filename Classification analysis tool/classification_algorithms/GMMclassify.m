function [labels, scores] = GMMclassify(gmm, data)

% [label, scores] = GMMclassify(gmm, data)
% 
% This function performs the classification of input data according to the
% Gaussian Mixture Models classifier.
% 
% gmm(i).class  - class ID
% gmm(i).obj    - object with the GMM model
% data          - input data matrix; columns are related to the features
%                 rows to the observations
% 
% label         - column vector with predicted numeric labels (classes)
% scores        - matrix containing the scores
%                 (each column is realted to one class

%% Paths and variables
scores = zeros(size(data, 1), size(gmm, 2));
labels = zeros(size(data, 1), 1);

%% Classify
for i = 1:size(data,1)
    for j = 1:size(gmm,2)
        scores(i, j) = pdf(gmm(j).obj, data(i, :));
    end
    
    [~, ind]  = max(scores(i, :));
    labels(i) = gmm(ind).class;
end