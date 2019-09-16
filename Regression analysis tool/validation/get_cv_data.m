function [data_sorted, labels_sorted, num_cat_vectors, cv_fold_sizes, ...
    random_order] = get_cv_data(data, labels, num_folds, categories)

% [data_sorted, labels_sorted, num_cat_vectors, cv_fold_sizes, ...
%    random_order] = get_cv_data(data, labels, num_folds, categories)
% 
% This function randomly prepares data for the cross-validation process
%
% data            - input data matrix (observations/features)
% labels          - class label stored in the associated input matrix
% categories      - column vector of class values in the input matrix
% num_folds       - number of cross-validation folds to preform
%
% data_sorted     - observations from data matrix grouped by category
%                   and in a random order within their category
% labels_sorted   - rows from categories vector grouped by category
%                   and in a random order within their category
% num_cat_vectors - column vector with the number of vectors in the data 
%                   set for each category (class)
% cv_fold_sizes   - matrix (with size [num_cat_vectors x num_folds] that
%                   contains the number of vectors to include in each of 
%                   the 'num_folds' folds for each category
% random_order    - random order (permutation of indices)
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

%% Get the number of vectors belonging to the categories (classes)
num_cat_vectors = get_num_category_vectors(labels, categories);

%% Get the cross-validation fold sizes for the provided data set
cv_fold_sizes = get_cv_fold_sizes(num_cat_vectors, num_folds);

%% Re-group the vectors according to categories (classes)
[data_sorted, labels_sorted, random_order] = ...
    rnd_sort_and_group(data, labels, categories);
