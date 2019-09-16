function [p_idx] = create_perm_idx(num_obs, num_perm, options)

% [p_idx] = create_perm_idx(num_obs, num_perm, sett)
%
% This function creates indices for the permutation test that is commonly
% used technique for the classification significance evaluation. Computed 
% indices are created randomly in range from [0, num_perm] and p_idx € N. 
%
% INPUTS:
% num_perm          - number of permutations (max. value of p_idx)
% options           - additional settings
%
% OUTPUTS:
% p_idx             - matrix of computed indices (rows: observations; 
%                     columns: permutations)
%
% SETTINGS STRUCTURE:
%   - options.save_to_mat       - switch to store the p_idx in *.mat file
%   - options.out_name          - output file name specification (*.mat)
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
if ((nargin < 3) || (isempty(options)))
    options.save_to_mat = true;  
    options.out_name    = [num2str(randi(intmax, 1, 1)) '.mat'];   
else
    if (~isfield(options, 'save_to_mat'))
        options.save_to_mat = true;
    end
    if (~isfield(options, 'out_name'))
        options.out_name = [num2str(randi(intmax, 1, 1)) '.mat'];  
    end
end

%% Create the random permutation (idx) matrix
p_idx = zeros(num_obs, num_perm);
for i = 1:num_perm
    p_idx(:, i) = randperm(num_obs).';
end

%% Store the data into *.mat files
if (options.save_to_mat)
    save(options.out_name, 'p_idx');
end