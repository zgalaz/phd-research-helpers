function [data_train, labels_train, data_test, labels_test] = ...
    get_fold_vectors(data_sorted, labels_sorted, categories,  ...
        num_cat_vectors, cv_fold_sizes, round_number)

% [data_train, labels_train, data_test, labels_test] =        ...
%    get_fold_vectors(data_sorted, labels_sorted, categories, ...
%        num_cat_vectors, cv_fold_sizes, round_number)
% 
% Selects the vectors to use for training and validation for the specified
% round number. With k-fold cross-validation, the data set is divided into
% k folds, and then training and validation are performed in 'k' separate
% rounds. In each round, 1 fold is used for validation while the remaining
% k - 1 folds are used for training.
%
% data_sorted     - input data matrix: (observations/features) grouped by 
%                   category values (classes), and in random order within
%                   their category -> out: rnd_sort_and_group()
% labels_sorted   - rows from categories vector grouped by category and in
%                   random order within their category
% categories      - column vector of class values in the input data matrix
% num_cat_vectors - column vector containing the number of vectors in data
%                   matrix belonging to each category (class)                  
% cv_fold_sizes   - matrix (with size [num_cat_vectors x num_folds] that
%                   contains the number of vectors to include in each of 
%                   the 'num_folds' folds for each category.
% round_number    - current round of cross-validation
%
% data_train      - training samples (sub-matrix for this cv round)
% labels_train    - training labels (sub-vector for this cv round)
% data_test       - testing samples (sub-matrix for this cv round)
% labels_test     - testing labels (sub-vector for this cv round) 
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

%% Create outpu vectors (empty at the begining)
data_train   = [];
data_test    = [];
labels_train = [];
labels_test  = [];

%% Verify that the vectors are properly sorted
cat_start = 1;

for cat = 1:length(categories)
    
    % Compute the index of the last vector of this category
    cat_end = cat_start + num_cat_vectors(cat) - 1;

    % Verify that all vectors in the range: <catStart:catEnd> 
    if (any(labels_sorted(cat_start:cat_end) ~= categories(cat)))
        error('Input vectors are not properly sorted');
    end
    
    % Set the starting index of the next category
    cat_start = cat_end + 1;
end

%% Get the number of folds from the foldSizes matrix
num_folds = size(cv_fold_sizes, 2);

%% Get the cross-validation fold vectors
cat_start = 1;

for cat = 1:length(categories)

    % Get the list of fold sizes for this category as a column vector
    cat_fold_sizes = cv_fold_sizes(cat, :)';
    
    % Set the starting index of the first cv fold for this category
    fold_start = cat_start;
    
    % Iterate over cross-validation folds and obtain the vectors
    for fold = 1:num_folds
        
        % Compute the index of the last vector in this fold
        fold_end = fold_start + cat_fold_sizes(fold) - 1;
        
        % Select all vectors in this fold
        fold_vectors = data_sorted(fold_start:fold_end, :);
        fold_classes = labels_sorted(fold_start:fold_end, :);
        
        % Append the vectors to the testing (validation) set
        if (fold == round_number)
            data_test   = [data_test;   fold_vectors];
            labels_test = [labels_test; fold_classes];
        else
            data_train   = [data_train;   fold_vectors];
            labels_train = [labels_train; fold_classes];
        end
        
        % Update the starting index of the next fold
        fold_start = fold_end + 1;
    end
    
    % Set the starting index of the next category
    cat_start = cat_start + num_cat_vectors(cat);   
end