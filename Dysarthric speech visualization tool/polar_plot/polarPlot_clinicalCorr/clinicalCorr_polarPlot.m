function out = clinicalCorr_polarPlot(features, clinScales, options)

% out = clinicalCorr_polarPlot(features, clinScales, options)
% 
% The function plots the correlation between the speech feature values and 
% the associated subjective clinical rating scales (UPDRS III, NMSS, etc.)
% of the patients with Parkinson's disease in a polar graph. The features 
% and also the rating scales can be defined as 3D matrices. The third dim.
% does represent measurements in time, which can be used to track changes
% in time over the measurements. 
%
% This function is designed to support multiple type of correlation coeff. 
% such as - Pearson's corr. coeff., Spearman's corr. coeff., Kendall's tau 
% and Goodman-Kruskal's gamma. The function does create a separate subplot 
% for each corr. coeff. to make the visualizations easier to read. It also 
% does check data integrity and performs the cleanup (removes observations
% with NaN, inf or complex values).
%
% Multiple graphical and computional settings can be changed according to
% user's preferences by manipulating the options structure. all the stuff
% that can be adjusted is listed bellow in options structure description.
%
% -------------------------------------------------------------------------
% TODO TASK LIST: 
% [1] Add visualization of p-values
% [2] Adjust the visualization of negative correlation
%     (now it displays the absolute values of computed corr. coeff.)
% 
% -------------------------------------------------------------------------
% INPUT VARIABLES:
% features             - 3D matrix of features 
%                           1. dimension: observations
%                           2. dimension: features
%                           3. dimension: measurements in time
% clinScales           - 3D matrix of clinical scales
%                           1. dimension: observations
%                           2. dimension: clinical scales
%                           3. dimension: measurements in time
% options              - settings (options) of graphical representation
%
% -------------------------------------------------------------------------
% OUTPUT VARIABLES:
% out                  - output structure that holds output information 
%                        stored in the following member variables:
%                           (timeInstances).Fh(feature).r - corr. coeff.
%                           (timeInstances).Fh(feature).p - p-values
%
% -------------------------------------------------------------------------
% OPTIONS STRUCTURE:
% options.typeCorrCoeff       - type of correlation coeff. to compute
% options.areaTransparency    - area transparency (alfa) parameter
% options.printCorr           - switch to display relative areas
% options.circleColorRGB      - color of the circles in the graph (RGB)
% options.lineColorRGB        - color of the lines in the graph (RGB!
% options.labelsFontType      - type of the font set to labels
% options.labelsFontWeight    - weight of the font set to labels
% options.labelsFontAngle     - angle of the font set to labels
% options.labelsFontSmooth    - smoothing of the font set to labels
% options.labelsFontSize      - size of the font set to labels
% options.labelsRotateDeg     - rotation of the font set to labels
% options.dataLineWidth       - width of the font set to data
% options.dataLineStyle       - style of the font set to data
% options.dataColorStyle      - color of the font set to data
% options.showLegend          - switch to show the legend (true/false)
% options.legendCellArray     - cell array defining the legend
% options.legendLocation      - location of the legend
% options.legendOrientation   - orientation of the legend
% options.legendBoxOutline    - switch to plot outline for legend
% options.legendEdgeColorRGB  - color of the edge of the legend
% options.legendTextColorRGB  - color of the text of the legend
% options.legendBGColorRGB    - color of the area of the legend
% options.legendLineWidth     - width of the line of the legend
% options.legendFontType      - type of the font of the legend
% options.legendFontWeight    - weight of the font of the legend
% options.legendFontAngle     - angle of the font of the legend
% options.legendFontSize      - size of the font of the legend
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

%% Pre-process input prosodic data
clinLabels = clinScales(1, :, 1);
clinData   = cell2mat(clinScales(2:end, :, 1));

if (size(features, 1) ~= size(clinData, 1))
    error('Observation does not match clinical data');
end

numObs    = size(features, 1);
numFeat   = size(features, 2);
numInst   = size(features, 3);
numScales = length(clinLabels);

%% Paths and variables
if ((nargin < 3) || (isempty(options)))
    options.areaTransparency   = 0.3;
    options.printCorr          = true;
    options.circleColorRGB     = [214/256 214/256 214/256];
    options.lineColorRGB       = [214/256 214/256 214/256];
    options.labelsFontType     = 'Times New Roman';
    options.labelsFontWeight   = 'bold';
    options.labelsFontAngle    = 'normal';
    options.labelsFontSmooth   = 'on';
    options.labelsFontSize     = 8;
    options.labelsRotateDeg    = 30;
    options.dataLineWidth      = 1.0;
    options.dataLineStyle      = '.-';
    options.dataColorStyle     = 'hsv';
    options.showLegend         = false;
    options.legendCellArray    = cellstr(num2str((1:numObs)'));
    options.legendLocation     = 'northwestoutside';
    options.legendOrientation  = 'vertical';
    options.legendBoxOutline   = 'off';
    options.legendEdgeColorRGB = [0.15 0.15 0.15];
    options.legendTextColorRGB = [0 0 0];
    options.legendBGColorRGB   = [1 1 1];
    options.legendLineWidth    = 0.5;
    options.legendFontType     = 'Times New Roman';
    options.legendFontWeight   = 'normal';
    options.legendFontAngle    = 'normal';
    options.legendFontSize     = 8;
    options.typeCorrCoeff      = 'Spearman';
                                  % {'Spearman', ...
                                  % 'Pearson',  ...
                                  % 'Kendall',  ...
                                  % 'Goodman-Kruskal'};
else
    
    if (~isfield(options, 'typeCorrCoeff'))
        options.typeCorrCoeff = 'Spearman';
                                  % {'Spearman', ...
                                  % 'Pearson',  ...
                                  % 'Kendall',  ...
                                  % 'Goodman-Kruskal'};
    end
    if (~isfield(options, 'areaTransparency'))
        options.areaTransparency = 0.3;
    end
    if (~isfield(options, 'printCorr'))
        options.printCorr = true;
    end
    if (~isfield(options, 'circircColorRGB'))
        options.circleColorRGB = [214/256 214/256 214/256];  
    end
    if (~isfield(options, 'lineColorRGB'))
        options.lineColorRGB = [214/256 214/256 214/256];  
    end
    if (~isfield(options, 'labelsFontType'))
        options.labelsFontType = 'Times New Roman';
    end
    if (~isfield(options, 'labelsFontWeight'))
        options.labelsFontWeight = 'normal';
    end
    if (~isfield(options, 'labelsFontAngle '))
        options.labelsFontAngle = 'normal';
    end
    if (~isfield(options, 'labelsFontSmooth'))
        options.labelsFontSmooth = 'on';
    end
    if (~isfield(options, 'labelsFontSize'))
        options.labelsFontSize = 8;
    end
    if (~isfield(options, 'labelsRotateDeg'))
        options.labelsRotateDeg = 30;
    end
    if (~isfield(options, 'dataLineWidth'))
        options.dataLineWidth = 1.0;
    end
    if (~isfield(options, 'dataLineStyle'))
        options.dataLineStyle = '.-';
    end
    if (~isfield(options, 'dataColorStyle'))
        options.dataColorStyle = 'hsv';
    end
    if (~isfield(options, 'showLegend'))
        options.showLegend = false;
    end
    if (~isfield(options, 'legendCellArray'))
        options.legendCellArray = cellstr(num2str((1:numObs)'));
    end
    if (~isfield(options, 'legendLocation'))
        options.legendLocation = 'northwestoutside';
    end
    if (~isfield(options, 'legendOrientation'))
        options.legendOrientation = 'vertical';
    end
    if (~isfield(options, 'legendBoxOutline'))
        options.legendBoxOutline = 'off';
    end
    if (~isfield(options, 'legendEdgeColorRGB'))
        options.legendEdgeColorRGB = [0.15 0.15 0.15];
    end
    if (~isfield(options, 'legendTextColorRGB'))
        options.legendTextColorRGB = [0 0 0];
    end
    if (~isfield(options, 'legendBGColorRGB'))
        options.legendBGColorRGB = [1 1 1];
    end
    if (~isfield(options, 'legendLineWidth'))
        options.legendLineWidth = 0.5;
    end
    if (~isfield(options, 'legendFontType'))
        options.legendFontType = 'Times New Roman';
    end
    if (~isfield(options, 'legendFontWeight'))
        options.legendFontWeight = 'normal';
    end
    if (~isfield(options, 'legendFontAngle'))
        options.legendFontAngle = 'normal';
    end
    if (~isfield(options, 'legendFontSize'))
        options.legendFontSize = 8;
    end
end

%% Prepare temporary data
corrSubPlot = false;
subplotMask = [1, 1, 1;   ...
               1, 2, 1;   ...
               2, 2, 1;   ...
               2, 2, 1;   ...
               3, 2, 1;   ...
               3, 2, 1;   ...
               3, 3, 1;   ...
               3, 3, 1;   ...
               3, 3, 3];
               
if (iscell(options.typeCorrCoeff))
    corrSubPlot  = true;
    numCorrCoeff = length(options.typeCorrCoeff);
    corrNumPlots = subplotMask(numCorrCoeff, :);
else
    numCorrCoeff = 1;
    corrNumPlots = subplotMask(1, :);
end

out = struct();

%% Process the data in time (if size(features, 3) > 1)
for time = 1:numInst
    data = features(:, :, time);
    clin = clinData(:, :, time);
    area = zeros(numFeat, numScales, numCorrCoeff);
    
    suplotCounter = 1;
    figure(time);
    
    % ---------------------------------------------------------------------
    % [A] Clean the data before processing
    %     Remove observation with NaN, inf or complex values to get only
    %     data that will not produce obscure values during a correlation
    sdel = false(numObs, 1);
    fdel = false(numObs, 1);
    
    for scale = 1:numScales
        sdel  = sdel | ((isnan(clin(:, scale)) |        ...
                isinf(clin(:, scale)) |                 ...
                imag(clin(:, scale))) == 1);
    end
    
    for feat = 1:numFeat
        fdel = fdel | ((isnan(data(:, feat)) |          ...
               isinf(data(:, feat)) |                   ...
               imag(data(:, feat))) == 1);
    end
    
    clin = real(clin(~(sdel | fdel), :));
    data = real(features(~(sdel | fdel), :));
    
    % ---------------------------------------------------------------------
    % [B] Perform the correlation ([clinical - paraclinical] data)
    %     Calculate the corr. coeff. and p-values for all chosen methods
    %     between the values of speech features, and associated clinical
    %     scales for iterated time instance (measurement).
    corr_r = zeros(numFeat, numScales, numCorrCoeff);
    corr_p = zeros(numFeat, numScales, numCorrCoeff);
    
    for coeff = 1:numCorrCoeff
        if (corrSubPlot)
            coeffType = options.typeCorrCoeff{coeff};
        else
            coeffType = options.typeCorrCoeff;
        end
 
        for scale = 1:numScales
            for feat = 1:numFeat
                [r, p] = perf_correlation(data(:, feat), ...
                    clin(:, scale),                      ...
                    coeffType);
                
                corr_r(feat, scale, coeff) = cell2mat(r);
                corr_p(feat, scale, coeff) = cell2mat(p);
            end
        end
    end
    
    % ---------------------------------------------------------------------
    % [C] Process the data for each correlation coefficient
    %     Plot the correlation graph representing the association between
    %     the feature values, and the associated clinical scales for each
    %     correlation coeff. into its own subplot.
    for coeff = 1:numCorrCoeff
        
        if (corrSubPlot)
            coeffType = options.typeCorrCoeff{coeff};
        else
            coeffType = options.typeCorrCoeff;
        end
        
        % -----------------------------------------------------------------
        % [01] Set the subplot
        %      Set the subplot for iterated correlation coefficient. Each
        %      coefficient will be graphed in its own subplot (for better
        %      readability of the plot).
        corrNumPlots(end) = suplotCounter;
        subplot(corrNumPlots(1), corrNumPlots(2), corrNumPlots(3));
        hold on;
        
        if (corrSubPlot)
            suplotCounter = +suplotCounter + 1;
        end
        
        % -----------------------------------------------------------------
        % [02] Set the circle parameters
        %      Set the basic set of parameters for all graphs, for example 
        %      radius and other position properties (center, etc.). Radius 
        %      is computed as the corr. coeff. achieved.
        radius = abs([corr_r(:, :, coeff) corr_r(:, 1, coeff)]');
        center = [0 0];
        circle = 2*pi;
        approx = 1000;

        % -----------------------------------------------------------------
        % [03] Draw outer circle
        %      Draw the outer circle that represents the maximum radius in 
        %      the circle computed for all correlation coefficients stored
        %      in the graph.
        amp = max(abs(radius(:)));
        phi = linspace(0, circle, approx)';

        x = amp*cos(phi) + center(1);
        y = amp*sin(phi) + center(2);
        plot(x, smooth(y), 'Color', options.circleColorRGB);
        
        
        % -----------------------------------------------------------------
        % [04] Draw middle-circles
        %      Draw the predefined number of middle-circles to make the 
        %      polar  plot easier to read. The number of middle-circles 
        %      is set into a steady variable: nCircles = 4 (this can be 
        %      adjusted if needed).
        nCircles = 4;
        kRadius  = linspace(0, amp, nCircles + 2);

        for rad  = 1:nCircles
            kAmp = kRadius(rad + 1);
            phi  = linspace(0, circle, approx)';

            x = kAmp*cos(phi) + center(1);
            y = kAmp*sin(phi) + center(2);
            plot(x, smooth(y), 'Color', options.circleColorRGB);
        end
        
        % -----------------------------------------------------------------
        % [05] Draw the actual series
        %      Draw the actual data in the graph. Connect the correlation  
        %      coefficients into the continuous graph differentiated with
        %      a different color (explanation is stored in the legend)
        phi = linspace(0, circle, numScales + 1)';
        x   = bsxfun(@times, radius, cos(phi) + center(1))';
        y   = bsxfun(@times, radius, sin(phi) + center(2))';
        mdl = zeros(1, numFeat);

        if (strcmpi(options.dataColorStyle, 'hsv'))
            clr = hsv(numFeat);
        end

        for dat = 1:numFeat
            mdl(dat) = plot(x(dat, :), y(dat, :),           ...
                options.dataLineStyle,                      ...
                'Color', clr(dat, :),                       ...
                'LineWidth', options.dataLineWidth);
        end
        
        % -----------------------------------------------------------------
        % [06] Draw lines (basic)
        %      Draw the default boundary lines that divide a whole circle 
        %      represented by 360 degrees into numFeat parts. These lines 
        %      are drawn bellow the graph and are set to gray color.
        phi = linspace(0, circle, numScales + 1)'; phi(end) = [];

        for lab = 1:numScales
            T_x = [0; amp]*cos(phi(lab)) - [0; 0]*sin(phi(lab));
            T_y = [0; amp]*sin(phi(lab)) + [0; 0]*cos(phi(lab));

            plot(T_x, T_y, 'Color', options.lineColorRGB);
        end
        
        % -----------------------------------------------------------------
        % [07] Compute the correlation areas
        %      Compute the absolute areas of the triangles that are formed  
        %      bycorrelations. It stores the information about the area of 
        %      each feature.
        for dat = 1:numScales
            for feat = 1:numFeat
                Xval = [0, x(feat, dat), x(feat, dat + 1)];
                Yval = [0, y(feat, dat), y(feat, dat + 1)];

                area(feat, dat, coeff) = +area(feat, dat, coeff) + ...
                    polyarea(Xval, Yval);
                
                h = fill(Xval, Yval, clr(feat, :));
                set(h, 'facealpha', options.areaTransparency);
            end
        end
        
        % -----------------------------------------------------------------
        % [08] Print the legend
        %      If the legend is enabled, the function displays the legend 
        %      in the predefined area with the settings stored inside the 
        %      'options' structure (e.g. legend text, colors, etc.).
        if (options.showLegend)
            legend(mdl, options.legendCellArray,            ...
                'Location', options.legendLocation,         ...
                'Orientation', options.legendOrientation,   ...
                'Box', options.legendBoxOutline,            ...
                'EdgeColor', options.legendEdgeColorRGB,    ...
                'TextColor', options.legendTextColorRGB,    ...
                'Color', options.legendBGColorRGB,          ...
                'LineWidth', options.legendLineWidth,       ...
                'FontName', options.legendFontType,         ...
                'FontWeight', options.legendFontWeight,     ...
                'FontAngle', options.legendFontAngle,       ...
                'FontSize', options.legendFontSize);
        end
        
        % -----------------------------------------------------------------
        % [09] Draw labels
        %      Draw the labels in the figure. Labels hold the information 
        %      about the clinical scales in the graph. Then function also 
        %      rotates the labels (if set to do so).
        phi = linspace(0, circle, numScales + 1)'; phi(end) = [];
        tmp = amp + 0.080;

        x = tmp*cos(phi) + center(1);
        y = tmp*sin(phi) + center(2);

        for lab = 1:numScales
            text(x(lab), y(lab), clinLabels{lab},           ...
                'FontName', options.labelsFontType,         ...
                'FontWeight', options.labelsFontWeight,     ...
                'FontSize', options.labelsFontSize,         ...
                'FontAngle', options.labelsFontAngle,       ...
                'FontSmoothing', options.labelsFontSmooth);
        end

        set(findobj(gca, 'Type', 'text'), ...
            'Rotation', options.labelsRotateDeg);
        
        % -----------------------------------------------------------------
        % [10] Set the axis properties
        %      Set the axis properties to fit the square sizes. Otherwise 
        %      the plot is non-symetric and needs to be manualy adjusted. 
        hold off;
        axis square;
        axis([-1.0 1.0 -1.0 1.0]*(amp + 0.2));
        title([num2str(time) '. time instance, '  ...
            'corr.: ' coeffType],                 ...
            'FontSize', options.labelsFontSize,   ...
            'FontName', options.labelsFontType);
    
        % -----------------------------------------------------------------
        % [11] Set the output variables
        %      Set the evolution of iterated correlation coefficients in 
        %      time for all the associated feature vectors, and clinical
        %      scales in the graph
        out(time).corr(coeff).r = corr_r;
        out(time).corr(coeff).p = corr_p;
    end
end