function num_cat_vectors = get_num_category_vectors(labels, categories)

% num_cat_vectors = get_num_category_vectors(labels, categories)
%
% This function returns the number of vectors stored in the data matrix
% belonging to each category. Categories that corresponds to the matrix 
% where rows correspond to the observations, and columns correnspond to
% features. Categories (classes) are stored in column vector 'labels'.
% The vectors may be supplied in any order (don't need to be grouped by
% category) and the category values do not have to be contiguous. 
%
% labels          - class label stored in the associated data matrix
% categories      - column vector of class values in the data matrix
%
% num_cat_vectors - column vector containing the number of vectors in
%                   data matrix: [observations/features] belonging to
%                   each category.
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

%% Get the number of categories (classes)
num_categories  = length(categories);

%% Create the vector to hold the number of categories (classes)
num_cat_vectors = zeros(num_categories, 1);

%% Check the format of classes (cell of strings/numeric array)
if (iscell(labels) && iscell(categories))
    labels     = conv_labels2mat(labels);
    categories = conv_labels2mat(categories);
end

%% Get the number of vectors belonging to the categories (classes)
for cat = 1:num_categories
    ctg = categories(cat);
    
    num_cat_vectors(cat, 1) = sum(labels == ctg);  
end