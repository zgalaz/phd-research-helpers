function [data_sorted, labels_sorted, random_order] = ...
    rnd_sort_and_group(data, labels, categories)

% [data_sorted, labels_sorted, random_order] = ...
%    rnd_sort_and_group(data, labels, categories)
% 
% This function randomly sort the vectors (observations) inside the data
% matrix and groups them by category. The function does the same process
% with labels (associated classes)
%
% data              - input data matrix (observations/features)
% labels            - class label stored in the associated input matrix
% categories        - column vector of class values in the input matrix
%
% data_sorted       - observations from data matrix grouped by category
%                     and in a random order within their category
% labels_sorted     - rows from categories vector grouped by category
%                     and in a random order within their category
% random_order      - random order (permutation of indices)
%
% Implemented according to:
% https://chrisjmccormick.wordpress.com/2013/07/31/ ... 
%   k-fold-cross-validation-with-matlab-code/
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

%% Get the total number of input vectors (observations)
total_vectors = size(data, 1);

%% Get a random order of the indeces
random_order  = randperm(total_vectors).';

%% Sort the vectors and categories with the random order
random_vectors = data(random_order, :);
random_classes = labels(random_order, :);

%% Create the empty vectors for sorted output variables
data_sorted   = [];
labels_sorted = [];

%% Re-group the vectors according to categories (classes)
for cat = 1:length(categories)
    ctg = categories(cat);
    
    % Select all of the vectors for this category
    cat_vectors = random_vectors((random_classes == ctg), :);
    cat_classes = random_classes((random_classes == ctg), :);
   
    % Append the vectors for this category
    data_sorted   = [data_sorted;   cat_vectors];
    labels_sorted = [labels_sorted; cat_classes];
end