function [out_results] = partial_correlation(feature_matrix,  ...
                                             clinical_matrix, ...
                                             covariates,      ...
                                             options)

% [out_results] = partial_correlation(feature_matrix,  ...
%                                     clinical_matrix, ...
%                                     covariates,      ...
%                                     options)
%
% This function computes the partial correlation between the data stored 
% in the input feature matrix and data stored in the input clinical data
% matrix, while controlling for the effect of covariates stored in the
% input covariates matrix.
%
% INPUTS:
% feature_matrix    - input feature matrix (numerical/cell array)
% clinical_matrix   - input clinical data (numerical/cell array)
% covariates        - input covariates data (numerical/cell array)
% options           - input settings (options)
%
% OUTPUTS:
% out_results       - output structure with the results
%                       .corr_results - results of the correlations
%                       .corr_tables  - results stored in tables
%						.lin_models   - generalized linear models
%
% OPTIONS:
%   - options.type_corr    - type/s of correlation coefficient/s
%   - options.analyse_BFBC - switch to perform BFBC* analysis
%   - options.analyse_BFDC - switch to perform BFDC** analysis
%   - options.analyse_DFDC - switch to perform DFDC*** analysis
%   - options.feat_labels  - labels for the features
%   - options.clin_labels  - labels for the clinical data
%   - options.covs_labels  - labels for the covariates
%   - options.create_xlsx  - switch to create *.xlsx output table/s
%
% *baseline features, baseline clinical data correlation
% **baseline features, delta clinical data correlation
% ***delta features, delta clinical data correlation
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

%% Check/process the input variables (mandatory)

% Check/process data in the feature matrix
if (iscell(feature_matrix))
    temp_matrix = feature_matrix;
    feat_matrix = zeros(            ... % Convert to numerical matrix:
        size(temp_matrix{1, 1}, 1), ... % 1. dim.: observations 
        size(temp_matrix{1, 1}, 2), ... % 2. dim.: features
        length(temp_matrix)         ... % 3. dim.: measurements (time)
    );
    
    % Create 3D feature matrix
    for i = 1:length(temp_matrix)
        feat_matrix(:, :, i) = real(temp_matrix{i});
    end
end

% Check/process data in the clinical matrix
if (iscell(clinical_matrix))
    temp_matrix = clinical_matrix;
    clin_matrix = zeros(            ... % Convert to numerical matrix:
        size(temp_matrix{1, 1}, 1), ... % 1. dim.: observations 
        size(temp_matrix{1, 1}, 2), ... % 2. dim.: features
        length(temp_matrix)         ... % 3. dim.: measurements (time)
    );
    
    % Create 3D clinical matrix
    for i = 1:length(temp_matrix)
        clin_matrix(:, :, i) = real(temp_matrix{i});
    end
end

% Check/process data in the covariates matrix
if (iscell(covariates))
    temp_matrix = covariates;
    covs_matrix = zeros(            ... % Convert to numerical matrix:
        size(temp_matrix{1, 1}, 1), ... % 1. dim.: observations 
        size(temp_matrix{1, 1}, 2), ... % 2. dim.: features
        length(temp_matrix)         ... % 3. dim.: measurements (time)
    );
    
    % Create 3D covariate matrix
    for i = 1:length(temp_matrix)
        covs_matrix(:, :, i) = real(temp_matrix{i});
    end
end

%% Check/process the input variables (optional)

% Check if the options are/not empty/filled
if ((nargin < 4) || isempty(options))
    options.type_corr    = {'Spearman', 'Pearson'};
    options.analyse_BFBC = true;
    options.analyse_BFDC = true;
    options.analyse_DFDC = true; 
    options.feat_labels  = cellstr(num2str((1:size(feat_matrix, 2))'));
    options.clin_labels  = cellstr(num2str((1:size(clin_matrix, 2))'));
    options.covs_labels  = cellstr(num2str((1:size(covs_matrix, 2))'));
    options.create_xlsx  = true;
else
    if (~isfield(options, 'type_corr'))
        options.type_corr = {'Spearman', 'Pearson'}; 
    end
    if (~isfield(options, 'analyse_BFBC'))
        options.analyse_BFBC = true;
    end
    if (~isfield(options, 'analyse_BFDC'))
        options.analyse_BFDC = true;
    end
    if (~isfield(options, 'analyse_DFDC'))
        options.analyse_DFDC = true;
    end
    if (~isfield(options, 'feat_labels'))
        options.feat_labels = cellstr(num2str((1:size(feat_matrix, 2))'));
    end
    if (~isfield(options, 'clin_labels'))
        options.clin_labels = cellstr(num2str((1:size(clin_matrix, 2))'));
    end
    if (~isfield(options, 'covs_labels'))
        options.covs_labels = cellstr(num2str((1:size(covs_matrix, 2))'));
    end
    if (~isfield(options, 'create_xlsx'))
        options.create_xlsx = true;
    end
end

%% Parse the input data properties

% Get the lengths/sizes of the input data
num_feat = size(feat_matrix, 2); % number of features
num_obs  = size(feat_matrix, 1); % number of observations
num_clin = size(clin_matrix, 2); % number of clinical features
num_covs = size(clin_matrix, 2); % number of covariates

% Check/process data in the options
if (iscell(options.type_corr))
    num_corr = length(options.type_corr);
    lab_corr = cell(1, num_corr*2);
    
    for i = 1:num_corr
        lab_corr{1, (i - 1)*num_corr + 1} = ...
            ['r (' options.type_corr{i} ')'];
        lab_corr{1, (i - 1)*num_corr + 2} = ...
            ['p (' options.type_corr{i} ')'];
    end
else
    num_corr = 1 ;
    lab_corr = cell(1, num_corr*2);
    
    lab_corr{1, 1} = ['r (' options.type_corr ')'];
    lab_corr{1, 2} = ['p (' options.type_corr ')'];
end

%% Correlate baseline features with baseline clinical data

% Prepare the output variables
corr_results_BFBC = cell(1, 2);
out_xlstable_BFBC = cell(1, 1);
out_linmodel_BFBC = cell(1, 1);

% Check if the BFBC analysis is required
if (options.analyse_BFBC)
    
    % Prepare the output variables
    corr_r  = zeros(num_feat, num_clin, num_corr);
    corr_p  = zeros(num_feat, num_clin, num_corr);
	lin_mdl = struct();
	
    % Iterate over all clinical data in clinical matrix
    for clin = 1:num_clin
        
        % Iterate over all features in feature matrix
        for feat = 1:num_feat
            
            % Parse the iterated data
            feat_data = feat_matrix(:, feat, 1);
            clin_data = clin_matrix(:, clin, 1);
            
            % Remove effect of covariates
            feat_data_mdl = fitlm(covs_matrix(:, :, 1), feat_data);
            clin_data_mdl = fitlm(covs_matrix(:, :, 1), clin_data);
            
			% Compute the residuals
            feat_data_residuals = feat_data_mdl.Residuals.Raw;
            clin_data_residuals = clin_data_mdl.Residuals.Raw;
            
			% Assign the linear models
			lin_mdl.feat_data_mdl(feat).lin_mdl = feat_data_mdl;
			lin_mdl.clin_data_mdl(clin).lin_mdl = clin_data_mdl;
	
            % Iterate over selected types of correlation
            for coeff = 1:num_corr
                func  = options.type_corr{coeff};
                
                % Compute the corr. coeff and p-value
                [r_coeff, p_coeff] = corr( ...
                    feat_data_residuals,   ...
                    clin_data_residuals,   ...
                    'type', func           ...
                );
                
                % Set the output variables
                corr_r(feat, clin, coeff) = r_coeff;
                corr_p(feat, clin, coeff) = p_coeff;
            end 
        end
    end
    
    % Merge the results
    corr_results_BFBC = {corr_r, corr_p};
    out_linmodel_BFBC = lin_mdl;
	
    % Store the results into *.xlsx table
    if (options.create_xlsx)
        
        % Iterate over all clinical data in clinical matrix
        for clin = 1:num_clin
            out_xlstable_BFBC = cell(num_feat + 1, num_corr*2 + 1);
            
			% Assign the labels
            out_xlstable_BFBC(1, :) = ...
                [{'features'}, lab_corr(:)'];
            out_xlstable_BFBC(:, 1) = ...
                [{'features'}; options.feat_labels(:)];
            
			% Compute the selected corr. coefficients
            for coeff = 1:num_corr
                data  = [corr_r(:, clin, coeff), corr_p(:, clin, coeff)];
                rows  = [2:size(out_xlstable_BFBC, 1)]';
                cols  = [((coeff - 1)*num_corr + 1 + 1) : ...
                    ((coeff - 1)*num_corr + 2 + 1)];
                
                out_xlstable_BFBC(rows, cols) = num2cell(data);
            end
            
            % Store the results into *.xlsx table
            table_name = 'partial_correlations_BFBC';
            sheet_name = options.clin_labels{clin};
            xlswrite([table_name '.xlsx'], out_xlstable_BFBC, sheet_name);
        end
    end
end

%% Correlate baseline features with delta clinical data

% Prepare the output variables
corr_results_BFDC = cell(1, 2);
out_xlstable_BFDC = cell(1, 1);
out_linmodel_BFDC = cell(1, 1);

% Check if the BFDC analysis is required
if ((options.analyse_BFDC) && ...
    (length(clinical_matrix) > 1))
    
    % Prepare the output variables
    corr_r  = zeros(num_feat, num_clin, num_corr);
    corr_p  = zeros(num_feat, num_clin, num_corr);
	lin_mdl = struct();

    % Iterate over all clinical data in clinical matrix
    for clin = 1:num_clin
        
        % Iterate over all features in feature matrix
        for feat = 1:num_feat
            
            % Parse the iterated data
            feat_data   = feat_matrix(:, feat, 1);
            clin_data_1 = clin_matrix(:, clin, 1);
            clin_data_2 = clin_matrix(:, clin, 2);
            
            % Remove effect of covariates
            feat_data_mdl   = fitlm(covs_matrix(:, :, 1), feat_data);
            clin_data_mdl_1 = fitlm(covs_matrix(:, :, 1), clin_data_1);
            clin_data_mdl_2 = fitlm(covs_matrix(:, :, 2), clin_data_2);
            
			% Compute the residuals
            feat_data_residuals   = feat_data_mdl.Residuals.Raw;
            clin_data_residuals_1 = clin_data_mdl_1.Residuals.Raw;
            clin_data_residuals_2 = clin_data_mdl_2.Residuals.Raw;
            delta_clin_residuals  = ...
                clin_data_residuals_2 - clin_data_residuals_1;
            
			% Assign the linear models
			lin_mdl.feat_data_mdl(feat).lin_mdl = feat_data_mdl;
			lin_mdl.clin_data_mdl(clin).lin_mdl = clin_data_mdl;
			
            % Iterate over selected types of correlation
            for coeff = 1:num_corr
                func  = options.type_corr{coeff};
                
                % Compute the corr. coeff and p-value
                [r_coeff, p_coeff] = corr( ...
                    feat_data_residuals,   ...
                    delta_clin_residuals,  ...
                    'type', func           ...
                );
                
                % Set the output variables
                corr_r(feat, clin, coeff) = r_coeff;
                corr_p(feat, clin, coeff) = p_coeff;
            end 
        end
    end
    
    % Merge the results
    corr_results_BFDC = {corr_r, corr_p};
    out_linmodel_BFDC = lin_mdl;
	
    % Store the results into *.xlsx table
    if (options.create_xlsx)
        
        % Iterate over all clinical data in clinical matrix
        for clin = 1:num_clin
            out_xlstable_BFDC = cell(num_feat + 1, num_corr*2 + 1);
            
			% Assign the labels
            out_xlstable_BFDC(1, :) = ...
                [{'features'}, lab_corr(:)'];
            out_xlstable_BFDC(:, 1) = ...
                [{'features'}; options.feat_labels(:)];
            
			% Compute the selected corr. coefficients
            for coeff = 1:num_corr
                data  = [corr_r(:, clin, coeff), corr_p(:, clin, coeff)];
                rows  = [2:size(out_xlstable_BFDC, 1)]';
                cols  = [((coeff - 1)*num_corr + 1 + 1) : ...
                    ((coeff - 1)*num_corr + 2 + 1)];
               
                out_xlstable_BFDC(rows, cols) = num2cell(data);
            end
            
            % Store the results into *.xlsx table
            table_name = 'partial_correlations_BFDC';
            sheet_name = options.clin_labels{clin};
            xlswrite([table_name '.xlsx'], out_xlstable_BFDC, sheet_name);
        end
    end
end

%% Correlate delta features with delta clinical data

% Prepare the output variables
corr_results_DFDC = cell(1, 2);
out_xlstable_DFDC = cell(1, 1);
out_linmodel_DFDC = cell(1, 1);

% Check if the DFDC analysis is required
if ((options.analyse_DFDC) && ...
    (length(feature_matrix) > 1) && ...
    (length(clinical_matrix) > 1))
    
    % Prepare the output variables
    corr_r  = zeros(num_feat, num_clin, num_corr);
    corr_p  = zeros(num_feat, num_clin, num_corr);
	lin_mdl = struct();
	
    % Iterate over all clinical data in clinical matrix
    for clin = 1:num_clin
        
        % Iterate over all features in feature matrix
        for feat = 1:num_feat
            
            % Parse the iterated data
            feat_data_1 = feat_matrix(:, feat, 1);
            feat_data_2 = feat_matrix(:, feat, 2);
            clin_data_1 = clin_matrix(:, clin, 1);
            clin_data_2 = clin_matrix(:, clin, 2);
            
            % Remove effect of covariates
            feat_data_mdl_1 = fitlm(covs_matrix(:, :, 1), feat_data_1);
            feat_data_mdl_2 = fitlm(covs_matrix(:, :, 2), feat_data_2);
            clin_data_mdl_1 = fitlm(covs_matrix(:, :, 1), clin_data_1);
            clin_data_mdl_2 = fitlm(covs_matrix(:, :, 2), clin_data_2);
            
			% Compute the residuals
            feat_data_residuals_1 = feat_data_mdl_1.Residuals.Raw;
            feat_data_residuals_2 = feat_data_mdl_2.Residuals.Raw;
            clin_data_residuals_1 = clin_data_mdl_1.Residuals.Raw;
            clin_data_residuals_2 = clin_data_mdl_2.Residuals.Raw;
            
            delta_feat_residuals  = ...
                feat_data_residuals_2 - feat_data_residuals_1;
            delta_clin_residuals  = ...
                clin_data_residuals_2 - clin_data_residuals_1;
            
			% Assign the linear models
			lin_mdl.feat_data_mdl(feat).lin_mdl = feat_data_mdl;
			lin_mdl.clin_data_mdl(clin).lin_mdl = clin_data_mdl;
			
            % Iterate over selected types of correlation
            for coeff = 1:num_corr
                func  = options.type_corr{coeff};
                
                % Compute the corr. coeff and p-value
                [r_coeff, p_coeff] = corr( ...
                    delta_feat_residuals,  ...
                    delta_clin_residuals,  ...
                    'type', func           ...
                );
                
                % Set the output variables
                corr_r(feat, clin, coeff) = r_coeff;
                corr_p(feat, clin, coeff) = p_coeff;
            end 
        end
    end
    
    % Merge the results
    corr_results_DFDC = {corr_r, corr_p};
    out_linmodel_DFDC = lin_mdl;
	
    % Store the results into *.xlsx table
    if (options.create_xlsx)
        
        % Iterate over all clinical data in clinical matrix
        for clin = 1:num_clin
            out_xlstable_DFDC = cell(num_feat + 1, num_corr*2 + 1);
            
			% Assign the labels
            out_xlstable_DFDC(1, :) = ...
                [{'features'}, lab_corr(:)'];
            out_xlstable_DFDC(:, 1) = ...
                [{'features'}; options.feat_labels(:)];
            
			% Compute the selected corr. coefficients
            for coeff = 1:num_corr
                data  = [corr_r(:, clin, coeff), corr_p(:, clin, coeff)];
                rows  = [2:size(out_xlstable_DFDC, 1)]';
                cols  = [((coeff - 1)*num_corr + 1 + 1) : ...
                    ((coeff - 1)*num_corr + 2 + 1)];
               
                out_xlstable_DFDC(rows, cols) = num2cell(data);
            end
            
            % Store the results into *.xlsx table
            table_name = 'partial_correlations_DFDC';
            sheet_name = options.clin_labels{clin};
            xlswrite([table_name '.xlsx'], out_xlstable_DFDC, sheet_name);
        end
    end
end

%% Merge/return the results
out_results = struct();

out_results.corr_results = [ ...
    corr_results_BFBC,       ...
    corr_results_BFDC,       ...
    corr_results_DFDC        ...
];

out_results.corr_tables = [ ...
    out_xlstable_BFBC,      ...
    out_xlstable_BFDC,      ...
    out_xlstable_DFDC       ...
];

out_results.lin_models = {  ...
	out_linmodel_BFBC,		...
	out_linmodel_BFDC,		...
	out_linmodel_DFDC		...
};