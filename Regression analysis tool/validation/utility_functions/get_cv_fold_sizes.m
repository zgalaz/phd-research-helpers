function cv_fold_sizes = get_cv_fold_sizes(num_cat_vectors, num_folds)

% cv_fold_sizes = get_cv_fold_sizes(num_cat_vectors, num_folds)
%
% Compute the cross-validation fold sizes for the provided data set. 
% num_cat_vectors - column vector with the number of vectors in the data 
%                   set for each category (class)
% num_folds       - number of folds used in the cross validation process
% cv_fold_sizes   - matrix (with size [num_cat_vectors x num_folds] that
%                   contains the number of vectors to include in each of 
%                   the 'num_folds' folds for each category
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
num_categories = length(num_cat_vectors);

%% For each category get the number of associated vectors
for cat = 1:num_categories
    
    % Get the number of associated vectors
    num_vectors = num_cat_vectors(cat, 1);
    
    % Verify that there are at least 'num_folds' samples
    if (num_vectors < num_folds)
        n = num2str(num_folds);
        error(['Each class must have at least ' n ' samples']);
    end
end

%% Create the output matrix of fold sizes
cv_fold_sizes = zeros(num_categories, num_folds);

%% Fill the output matrix (cv fold sizes)
% NOTE: cv_fold_sizes: matrix holding the number of vectors to place in 
%       each fold for each category. The number of folds may not divide
%       evenly into the number of vectors, so we need to distribute the
%       remainder.
for cat = 1:num_categories
    
    % Get the number of associated vectors
    num_vectors = num_cat_vectors(cat, 1);
            
    % Fill the output matrix iterating over all cv folds
    for fold = 1:num_folds
        
        % Divide the remaining number of vectors by the remaining 
        % number of folds (in order to distribute the remainder)
        % and store the fold size into the output matrix
        fold_size = ceil(num_vectors/(num_folds - fold + 1));
        cv_fold_sizes(cat, fold) = fold_size;
        
        % Update the number of remaining vectors for this category
        num_vectors = num_vectors - fold_size;
    end
end

%% Verify the fold sizes sum up correctly
if (any(sum(cv_fold_sizes, 2) ~= num_cat_vectors))
    error('Sum of fold sizes ~= the number of category vectors');
end