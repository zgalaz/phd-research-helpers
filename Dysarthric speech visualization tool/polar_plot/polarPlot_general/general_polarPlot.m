function out = general_polarPlot(data, labels, options)

% out = general_polarPlot(data, labels, options)
% 
% The function plots input speech feature features into a polar graph. The
% graph can plot the features for several speakers, and also for a several 
% measurements over time. Input data does contain the 3D matrix of speech 
% features and the associated labels and the options.
%
% -------------------------------------------------------------------------
% INPUT VARIABLES:
% data                 - structure with 3D matrices of speech features
%                           - 1.dim: observations (PD/HC etc.)
%                           - 2.dim: features; (e.g. std(F0), etc.)
%                           - 3.dim: measurements in time
% labels               - cell array of labels (for each observation)
% options              - settings (options) of graphical representation
%
% -------------------------------------------------------------------------
% OUTPUT VARIABLES:
% out                  - output structure that holds output information
%
% -------------------------------------------------------------------------
% OPTIONS STRUCTURE:
% options.normalizeData       - switch to turn [on/off] normalization
% options.normalizationMethod - method used to normalize data
% options.areaTransparency    - area transparency (alfa) parameter
% options.printArea           - switch to display relative areas
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
% ing. Zolt�n Gal�
% xgalaz00@stud.feec.vutbr.cz       
% 
% Department of Telecommunications
% Faculty of Electrical Engineering and Communication
% Brno University of Technology

%% Pre-process input prosodic data
numObs  = size(data, 1);
numFeat = size(data, 2);
numInst = size(data, 3);

%% Paths and variables
if ((nargin < 3) || (isempty(options)))
    options.normalizeData       = true;
    options.normalizationMethod = 'SoftMax';
    options.areaTransparency    = 0.3;
    options.printArea           = true;
    options.circleColorRGB      = [214/256 214/256 214/256];
    options.lineColorRGB        = [214/256 214/256 214/256];
    options.labelsFontType      = 'Times New Roman';
    options.labelsFontWeight    = 'bold';
    options.labelsFontAngle     = 'normal';
    options.labelsFontSmooth    = 'on';
    options.labelsFontSize      = 8;
    options.labelsRotateDeg     = 30;
    options.dataLineWidth       = 1.0;
    options.dataLineStyle       = '.-';
    options.dataColorStyle      = 'hsv';
    options.showLegend          = false;
    options.legendCellArray     = cellstr(num2str((1:numObs)'));
    options.legendLocation      = 'northwestoutside';
    options.legendOrientation   = 'vertical';
    options.legendBoxOutline    = 'off';
    options.legendEdgeColorRGB  = [0.15 0.15 0.15];
    options.legendTextColorRGB  = [0 0 0];
    options.legendBGColorRGB    = [1 1 1];
    options.legendLineWidth     = 0.5;
    options.legendFontType      = 'Times New Roman';
    options.legendFontWeight    = 'normal';
    options.legendFontAngle     = 'normal';
    options.legendFontSize      = 8;
else
    
    if (~isfield(options, 'normalizeData'))
        options.normalizeData = true; 
    end
    if (~isfield(options, 'normalizationMethod'))
        options.normalizationMethod = 'SoftMax'; 
    end
    if (~isfield(options, 'areaTransparency'))
        options.areaTransparency = 0.3;
    end
    if (~isfield(options, 'printArea'))
        options.printArea = true;
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

if ((nargin < 2) || (isempty(labels)))
    labels = cellstr(num2str((1:numFeat)'));
end

if (size(data, 2) ~= length(labels))
    error('features does not match the labels');
end

dataMtx = data;
out     = struct();

%% Process the data in time (if size(dataMtx, 3) > 1)
for time = 1:numInst
    data = dataMtx(:, :, time);
    
    figure(time);
    hold on;
    
    % ---------------------------------------------------------------------
    % [A] Set the circle parameters
    %     Set the basic set of parameters of the graph, such as its radius
    %     and other position properties (center, etc.). Radius is computed
    %     using a normalization method (if enabled).
    if (options.normalizeData)
        alg    = @normalize;
        radius = alg([data data(:, 1)]', options.normalizationMethod);
    else
        radius = [data data(:, 1)]';
    end
    
    center = [0 0];
    circle = 2*pi;
    approx = 1000;
    
    % ---------------------------------------------------------------------
    % [B] Draw outer circle
    %     Draw the outer circle denoting the maximum radius in the circle
    %     computed from all features stored in the graph. 
    amp = max(abs(radius(:)));
    phi = linspace(0, circle, approx)';
    
    x = amp*cos(phi) + center(1);
    y = amp*sin(phi) + center(2);
    plot(x, smooth(y), 'Color', options.circleColorRGB);

    % ---------------------------------------------------------------------
    % [C] Draw middle-circles
    %     Draw predefined number of middle-circles to make the polar plot
    %     easier to read. the number of middle plots is set to a steady
    %     variable: nCircles = 4 (can be adjusted if needed).
    nCircles = 4;
    kRadius  = linspace(0, amp, nCircles + 2);
    
    for rad  = 1:nCircles
        kAmp = kRadius(rad + 1);
        phi  = linspace(0, circle, approx)';

        x = kAmp*cos(phi) + center(1);
        y = kAmp*sin(phi) + center(2);
        plot(x, smooth(y), 'Color', options.circleColorRGB);
    end

    % ---------------------------------------------------------------------
    % [D] Draw the actual series
    %     Draw the actual data in the graph. Connect the features into the
    %     continuous graph. Data for each observation (speakers) is drawn
    %     in a different color (explanation is stored in the legend).
    phi = linspace(0, circle, numFeat + 1)';
    x   = bsxfun(@times, radius, cos(phi) + center(1))';
    y   = bsxfun(@times, radius, sin(phi) + center(2))';
    mdl = zeros(1, numObs);

    if (strcmpi(options.dataColorStyle, 'hsv'))
        clr = hsv(numObs);
    end

    for dat = 1:numObs
        mdl(dat) = plot(x(dat, :), y(dat, :),           ...
            options.dataLineStyle,                      ...
            'Color', clr(dat, :),                       ...
            'LineWidth', options.dataLineWidth);
    end
    
    % ---------------------------------------------------------------------
    % [E] Draw lines
    %     Draw the default boundary lines that divide the whole circle 360
    %     degrees into numFeat parts. The lines are drawn bellow the graph 
    %     and are set to gray color.
    phi  = linspace(0, circle, numFeat + 1)'; phi(end) = [];
    
    for lab = 1:numFeat
        T_x = [0; amp]*cos(phi(lab)) - [0; 0]*sin(phi(lab));
        T_y = [0; amp]*sin(phi(lab)) + [0; 0]*cos(phi(lab));
        
        plot(T_x, T_y, 'Color', options.lineColorRGB);
    end
    
    % ---------------------------------------------------------------------
    % [F] Compute the observation areas
    %     Calculate the absolute areas of the triangles formed by values
    %     of all features of the observatinos. It stores the information
    %     about the area of the features for all observations (speakers)
    area = zeros(numObs, numFeat);
    
    for dat = 1:numObs
        for feat = 1:numFeat
            Xval = [0, x(dat, feat), x(dat, feat + 1)];
            Yval = [0, y(dat, feat), y(dat, feat + 1)];
            
            area(dat, feat) = +area(dat, feat) + polyarea(Xval, Yval);
            h = fill(Xval, Yval, clr(dat, :));
            set(h, 'facealpha', options.areaTransparency);
        end
    end
    
    % ---------------------------------------------------------------------
    % [G] Print the legend
    %     If the legend's enabled, the function display the legend in the
    %     predefined area with the settings stored inside the 'options' 
    %     structure (e.g. legend text, colors, font settings, etc.).
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

    % ---------------------------------------------------------------------
    % [H] Draw labels
    %     Draw the labels in the figure. Labels does hold the information 
    %     about the features in the graph. This function also rotates the 
    %     labels (if set to do so).
    phi = linspace(0, circle, numFeat + 1)'; phi(end) = [];
    tmp = amp + 0.080;

    x = tmp*cos(phi) + center(1);
    y = tmp*sin(phi) + center(2);

    for lab = 1:numFeat
        text(x(lab), y(lab), labels{lab},               ...
            'FontName', options.labelsFontType,         ...
            'FontWeight', options.labelsFontWeight,     ...
            'FontSize', options.labelsFontSize,         ...
            'FontAngle', options.labelsFontAngle,       ...
            'FontSmoothing', options.labelsFontSmooth);
    end
    
    set(findobj(gca, 'Type', 'text'), 'Rotation', options.labelsRotateDeg);
        
    % ---------------------------------------------------------------------
    % [I] Set the axis properties
    %     Set the axis properties to fit the square sizes. Otherwise the 
    %     plot is non-symetric and it needs to be manualy adjusted. 
    hold off;
    axis square;
    axis([-1.0 1.0 -1.0 1.0]*(amp + 0.2));
    title([num2str(time) '. time instance'],  ...
        'FontSize', options.labelsFontSize,   ...
        'FontName', options.labelsFontType);
    
    % ---------------------------------------------------------------------
    % [J] Set the output variables
    %     Set the evolution of relative areas of all three prosodic models
    %     in time the associated feature vectors for all observations (the
    %     observations - speakers, e.g. (PD/HC)). 
    for dat = 1:numObs
        out(time).relArea(dat) = sum(area(dat, :))/(pi*amp^2); 
    end
end

%% Compute relative areas for feature groups
relArea = zeros(numObs, numInst);

for time = 1:numInst
    for dat = 1:numObs
        relArea(dat, time) = out(time).relArea(dat);
    end
end

%% Pllot the evolution of relative areas in time
if (numInst > 1)
    figure(time + 1)
    
    for dat = 1:numObs
        plot(relArea(dat, :),                           ...
            options.dataLineStyle,                      ...
            'Color', clr(dat, :),                       ...
            'LineWidth', options.dataLineWidth);
        
        xlabel('evolution of monopitch features',       ...
            'FontSize', options.labelsFontSize,         ...
            'FontName', options.labelsFontType);
        ylabel('feature values',                        ...
            'FontSize', options.labelsFontSize,         ...
            'FontName', options.labelsFontType);
        
        hold on;
        grid on; 
    end
end