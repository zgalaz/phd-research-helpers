function out = dysarthric_polarPlot(dysarthria, labels, options)

% out = dysarthric_polarPlot(dysarthria, labels, options)
% 
% The function plots dysarthric features (area of: phonation, articulation, 
% prosody and speech fluency) into a polar graph. The graph can display the 
% features for several speakers and also for several measurements over time
% Input data does contain the 3D matrices of dysarthric speech features:
%   -   phon (Ph - phonation)
%   -   artc (Ar - articulation)
%   -   pros (Pr - prosody)
%   -   rate (Sr - speech rate)
%
% -------------------------------------------------------------------------
% INPUT VARIABLES:
% prosody              - structure with 3D matrices for dysarthric features
%                        for Ph, Ar, Pr and Sr features, stored in:
%                           dysarthria.Ph_features
%                           dysarthria.Ar_features
%                           dysarthria.Pr_features
%                           dysarthria.Sr_features
%                        
%                        - 1.dim: observations (PD/HC etc.) 
%                        - 2.dim: features; (e.g. std(TEO), etc.)
%                        - 3.dim: measurements in time
%
% labels               - cell array of labels (for each feature)
% options              - settings (options) of graphical representation
%
% -------------------------------------------------------------------------
% OUTPUT VARIABLES:
% out                  - output structure that holds output information 
%                        stored in the following member variables:
%                           (timeInstances).Fh(observation).features
%                           (timeInstances).Ar(observation).features
%                           (timeInstances).Pr(observation).features
%                           (timeInstances).Sr(observation).features
%                           (timeInstances).Fh(observation).relArea
%                           (timeInstances).Ar(observation).relArea
%                           (timeInstances).Pr(observation).relArea
%                           (timeInstances).Sr(observation).relArea
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
% ing. Zoltán Galáž
% xgalaz00@stud.feec.vutbr.cz       
% 
% Department of Telecommunications
% Faculty of Electrical Engineering and Communication
% Brno University of Technology

%% Pre-process input dysarthric data
Ph = [];
Ar = [];
Pr = [];
Sr = [];

if (isfield(dysarthria, 'Ph_features'))
    Ph = dysarthria.Ph_features;
end
if (isfield(dysarthria, 'Ar_features'))
    Ar = dysarthria.Ar_features;
end
if (isfield(dysarthria, 'Pr_features'))
    Pr = dysarthria.Pr_features;
end
if (isfield(dysarthria, 'Sr_features'))
    Sr = dysarthria.Sr_features;
end

str = 'Dysarthric submodels';
s1  = [size(Ph, 1) size(Ar, 1) size(Pr, 1) size(Sr, 1)];
s2  = [size(Ph, 2) size(Ar, 2) size(Pr, 2) size(Sr, 2)];
s3  = [size(Ph, 3) size(Ar, 3) size(Pr, 3) size(Sr, 3)];
s1  = s1(s1 ~= 0);
s2  = s2(s2 ~= 0);
s3  = s3(s3 ~= 1);
L   = s2;

if (sum(s1 == s1(1)) ~= length(s1))
    error([str ' must have the same number of observations']);
end
if (sum(s3 == s3(1)) ~= length(s3))
    error([str ' must have the same number of instances']);
end

dataMtx = [Ph Ar Pr Sr];
numObs  = size(dataMtx, 1);
numFeat = size(dataMtx, 2);
numInst = size(dataMtx, 3);
numDys  = sum(~cellfun('isempty', struct2cell(dysarthria)));

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
    options.labelsFontSize      = 10;
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
    options.legendFontSize      = 10;
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
        options.labelsFontSize = 10;
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
        options.legendFontSize = 10;
    end
end

if ((nargin < 2) || (isempty(labels)))
    labels = cellstr(num2str((1:numFeat)'));
end

if (size(dataMtx, 2) ~= length(labels))
    error('features does not match the labels');
end

out = struct();

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
    
    center  = [0 0];
    circle  = 2*pi;
    partial = circle/numDys;
    approx  = 1000;
    
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
    phi = [];
    for dys = 1:numDys
        if (dys ~= numDys)
            tmp = linspace(partial*(dys - 1), partial*dys, L(dys) + 1)';
            phi = [phi; tmp(1: end - 1)];
        else
            tmp = linspace(partial*(dys - 1), partial*dys, L(dys) + 1)';
            phi = [phi; tmp];
        end
    end
    
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
    % [E] Draw lines (basic)
    %     Draw the default boundary lines that divide the whole circle 360
    %     degrees into numFeat parts. The lines are drawn bellow the graph 
    %     and are set to gray color.
    phi = [];
    for dys = 1:numDys
        if (dys ~= numDys)
            tmp = linspace(partial*(dys - 1), partial*dys, L(dys) + 1)';
            phi = [phi; tmp(1: end - 1)];
        else
            tmp = linspace(partial*(dys - 1), partial*dys, L(dys) + 1)';
            phi = [phi; tmp];
        end
    end
    
    phi(end) = [];
    plotBins = cumsum([1 L]);
    
    for lab = 1:numFeat
        T_x = [0; amp]*cos(phi(lab)) - [0; 0]*sin(phi(lab));
        T_y = [0; amp]*sin(phi(lab)) + [0; 0]*cos(phi(lab));
        
        if (~any(plotBins == lab))
            plot(T_x, T_y, 'Color', options.lineColorRGB);
        end
    end
    
    % ---------------------------------------------------------------------
    % [F] Compute the prosodic areas
    %     Compute the absolute areas of the triangles formed by prosodic
    %     features. It stores the information about the area of (Ph, Ar 
    %     Pr and Sr) features for all observations (speakers).
    area = zeros(numObs, numDys);
    bins = cumsum([1 L]);
    bins = bins(2:end);
    
    for dat  = 1:numObs
        iter = 1;
    
        for feat = 1:numFeat
            Xval = [0, x(dat, feat), x(dat, feat + 1)];
            Yval = [0, y(dat, feat), y(dat, feat + 1)];
            
            if (feat > bins(iter))
                iter = +iter + 1;
            end
            
            area(dat, iter) = +area(dat, iter) + polyarea(Xval, Yval);
            h = fill(Xval, Yval, clr(dat, :));
            set(h, 'facealpha', options.areaTransparency);
        end
    end
    
    % ---------------------------------------------------------------------
    % [G] Draw lines (prosodic borders)
    %     Draw the prosodic boundary lines to divide the whole circle 360
    %     degrees into three parts (Ph, Ar, Pr, Sr). These lines are drawn 
    %     on top of the graph and are bolder than default ones.
    phi = [];
    for dys = 1:numDys
        if (dys ~= numDys)
            tmp = linspace(partial*(dys - 1), partial*dys, L(dys) + 1)';
            phi = [phi; tmp(1: end - 1)];
        else
            tmp = linspace(partial*(dys - 1), partial*dys, L(dys) + 1)';
            phi = [phi; tmp];
        end
    end
    
    phi(end) = [];
    plotBins = cumsum([1 L]);
    
    for lab = 1:numFeat
        T_x = [0; amp]*cos(phi(lab)) - [0; 0]*sin(phi(lab));
        T_y = [0; amp]*sin(phi(lab)) + [0; 0]*cos(phi(lab));
        
        if (any(plotBins == lab))
            plot(T_x, T_y, '-k', 'LineWidth', 1.5);
        end
    end
     
    % ---------------------------------------------------------------------
    % [H] Print the legend
    %     If the legend is enabled, the function display the legend in the
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
    % [I] Draw labels
    %     Draw the labels in the figure. Labels hold the information about
    %     the prosodic features in the graph. The function also displays
    %     informative text explaining the feature areas (Ph, Ar, Pr, Sr) 
    %     and rotates the labels (if set to do so).
    phi = [];
    for dys = 1:numDys
        if (dys ~= numDys)
            tmp = linspace(partial*(dys - 1), partial*dys, L(dys) + 1)';
            phi = [phi; tmp(1: end - 1)];
        else
            tmp = linspace(partial*(dys - 1), partial*dys, L(dys) + 1)';
            phi = [phi; tmp];
        end
    end
    
    phi(end) = [];
    tmpAmp   = amp + 0.080;

    x = tmpAmp*cos(phi) + center(1);
    y = tmpAmp*sin(phi) + center(2);

    for lab = 1:numFeat
        text(x(lab), y(lab), labels{lab},               ...
            'FontName', options.labelsFontType,         ...
            'FontWeight', options.labelsFontWeight,     ...
            'FontSize', options.labelsFontSize,         ...
            'FontAngle', options.labelsFontAngle,       ...
            'FontSmoothing', options.labelsFontSmooth);
    end

    set(findobj(gca, 'Type', 'text'), 'Rotation', options.labelsRotateDeg);
    plotDys = true;
    cnt = 1;
    
    switch (numDys)
        case 4
            posX = [0.69 -0.97 -0.97  0.59];
            posY = [0.95  0.95 -0.95 -0.95]; 
        case 3
            posX = [-0.97 -0.97 -0.97];
            posY = [-0.75 -0.85 -0.95];      
        case 2
            posX = [-0.97 -0.97];
            posY = [ 0.95 -0.95]; 
        otherwise
            plotDys = false;
    end
    
    if (plotDys)
        if (~isempty(Ph))
            text(posX(cnt), posY(cnt), 'Phonation',         ...
                'FontSize', options.labelsFontSize,         ...
                'FontName', options.labelsFontType);
            cnt = +cnt + 1;
        end
        if (~isempty(Ar))
            text(posX(cnt), posY(cnt), 'Articulation',      ...
                'FontSize', options.labelsFontSize,         ...
                'FontName', options.labelsFontType);
            cnt = +cnt + 1;
        end
        if (~isempty(Pr))
            text(posX(cnt), posY(cnt), 'Prosody',           ...
                'FontSize', options.labelsFontSize,         ...
                'FontName', options.labelsFontType);
            cnt = +cnt + 1;
        end
        if (~isempty(Sr))
            text(posX(cnt), posY(cnt), 'Speech fluency',    ...
                'FontSize', options.labelsFontSize,         ...
                'FontName', options.labelsFontType);
        end
    end
         
    % ---------------------------------------------------------------------
    % [J] Set the axis properties
    %     Set the axis properties to fit the square sizes. Otherwise the
    %     plot is non-symetric and it needs to be manualy adjusted.
    hold off;
    axis square;
    axis([-1.0 1.0 -1.0 1.0]*(amp + 0.2));
    title([num2str(time) '. time instance'],  ...
        'FontSize', options.labelsFontSize,   ...
        'FontName', options.labelsFontType);
    
    % ---------------------------------------------------------------------
    % [K] Set the output variables
    %     Set the evolution of relative areas of all three prosodic models
    %     in time the associated feature vectors for all observations (the
    %     observations - speakers, e.g. (PD/HC)).
    for dat = 1:numObs
        cnt = 1;
        den = ((pi*amp^2)/numDys);
        
        if (~isempty(Ph))
            out(time).Ph(dat).relArea  = area(dat, cnt)/den; 
            out(time).Ph(dat).features = Ph(dat, :, time);
            cnt = +cnt + 1;
        end
        if (~isempty(Ar))
            out(time).Ar(dat).relArea  = area(dat, cnt)/den;
            out(time).Ar(dat).features = Ar(dat, :, time);
            cnt = +cnt + 1;
        end
        if (~isempty(Pr))
            out(time).Pr(dat).relArea  = area(dat, cnt)/den;
            out(time).Pr(dat).features = Pr(dat, :, time);
            cnt = +cnt + 1;
        end
        if (~isempty(Sr))
            out(time).Sr(dat).relArea = area(dat, cnt)/den;
            out(time).Sr(dat).features = Sr(dat, :, time);
        end
    end
end

%% Compute relative areas for feature groups
if (~isempty(Ph))
    relAreaPh = zeros(numObs, numInst);
end
if (~isempty(Ar))
    relAreaAr = zeros(numObs, numInst);
end
if (~isempty(Pr))
    relAreaPr = zeros(numObs, numInst);
end
if (~isempty(Sr))
    relAreaSr = zeros(numObs, numInst);
end

for time = 1:numInst
    for dat = 1:numObs
        if (~isempty(Ph))
            relAreaPh(dat, time) = out(time).Ph(dat).relArea;
        end
        if (~isempty(Ar))
            relAreaAr(dat, time) = out(time).Ar(dat).relArea;
        end
        if (~isempty(Pr))
            relAreaPr(dat, time) = out(time).Pr(dat).relArea;
        end
        if (~isempty(Sr))
            relAreaSr(dat, time) = out(time).Sr(dat).relArea;
        end
    end
end

%% Print relative areas
D = '-------------------------------------------------------------------';
    
if (options.printArea)   
    fprintf('\n'                                                        );
    fprintf('Relative dysarthric areas:\n'                              );
        
    for time = 1:numInst
        frm  = '';
        
        if (~isempty(Ph))
            frm = [frm 'Phonation' '\t'];
        end
        if (~isempty(Ar))
            frm = [frm 'Articulation' '\t'];
        end
        if (~isempty(Pr))
            frm = [frm 'Prosody' '\t' '\t'];
        end
        if (~isempty(Sr))
            frm = [frm 'Fluency' '\t' '\t'];
        end
        
        frm  = [frm 'Speaker'];
        fprintf('%s', D                                                 );
        fprintf('\n'                                                    );
        fprintf('%d. time instance (measurement):\n', time              );
        fprintf('\n'                                                    );
        fprintf(['\t' frm '\n']                                         );
        
        for dat = 1:numObs
            obs = options.legendCellArray{dat};
            
            if (~isempty(Ph))
                fprintf('\t%2.4f', relAreaPh(dat, time));
            end
            if (~isempty(Ar))
                fprintf('\t\t%2.4f', relAreaAr(dat, time));
            end
            if (~isempty(Pr))
                fprintf('\t\t\t%2.4f', relAreaPr(dat, time));
            end
            if (~isempty(Sr))
                fprintf('\t\t%2.4f', relAreaSr(dat, time));
            end
        
            fprintf('\t\t%s\n', obs); 
        end
    end 
    fprintf('\n');
end

%% Plot the evolution of relative areas in time
if (numInst > 1)
    figure(time + 1)
    plotSub = true;
    cnt = 1;
    
    switch (numDys)
        case 4
            subplotMtx = [2 2]; 
        case 3
            subplotMtx = [3 1];  
        case 2
            subplotMtx = [2 1];
        otherwise
            plotSub = false;
    end
    
    % ---------------------------------------------------------------------
    % [A] Plot relative area of Ph evolution in time
    if (~isempty(Ph))
        if (plotSub)
            subplot(subplotMtx(1), subplotMtx(2), cnt);
        end
        
        for dat = 1:numObs 
            plot(relAreaPh(dat, :),                             ...
                options.dataLineStyle,                          ...
                'Color', clr(dat, :),                           ...
                'LineWidth', options.dataLineWidth);

            xlabel('evolution of prosodic features',            ...
                'FontSize', options.labelsFontSize,             ...
                'FontName', options.labelsFontType);
            ylabel('feature values',                            ...
                'FontSize', options.labelsFontSize,             ...
                'FontName', options.labelsFontType);

            hold on;
            grid on; 
        end
        
        cnt = +cnt + 1;
    end
    
    % ---------------------------------------------------------------------
    % [B] Plot relative area of Ar evolution in time
    if (~isempty(Ar))
        if (plotSub)
            subplot(subplotMtx(1), subplotMtx(2), cnt);
        end
        
        for dat = 1:numObs
            plot(relAreaAr(dat, :),                             ...
                options.dataLineStyle,                          ...
                'Color', clr(dat, :),                           ...
                'LineWidth', options.dataLineWidth);

            xlabel('evolution of articulation features',        ...
                'FontSize', options.labelsFontSize,             ...
                'FontName', options.labelsFontType);
            ylabel('feature values',                            ...
                'FontSize', options.labelsFontSize,             ...
                'FontName', options.labelsFontType);

            hold on;
            grid on;
        end
        
        cnt = +cnt + 1;
    end
    
    % ---------------------------------------------------------------------
    % [C] Plot relative area of Pr evolution in time
    if (~isempty(Pr))
        if (plotSub)
            subplot(subplotMtx(1), subplotMtx(2), cnt);
        end
        
        for dat = 1:numObs
            plot(relAreaPr(dat, :),                         ...
                options.dataLineStyle,                      ...
                'Color', clr(dat, :),                       ...
                'LineWidth', options.dataLineWidth);

            xlabel('evolution of prosodic features',        ...
                'FontSize', options.labelsFontSize,         ...
                'FontName', options.labelsFontType);
            ylabel('feature values',                        ...
                'FontSize', options.labelsFontSize,         ...
                'FontName', options.labelsFontType);

            hold on;
            grid on;
        end
        cnt = +cnt + 1;
    end
    
    % ---------------------------------------------------------------------
    % [D] Plot relative area of Sr evolution in time
    if (~isempty(Sr))
        if (plotSub)
            subplot(subplotMtx(1), subplotMtx(2), cnt);
        end
        
        for dat = 1:numObs
            plot(relAreaSr(dat, :),                         ...
                options.dataLineStyle,                      ...
                'Color', clr(dat, :),                       ...
                'LineWidth', options.dataLineWidth);

            xlabel('evolution of speech fluency features',  ...
                'FontSize', options.labelsFontSize,         ...
                'FontName', options.labelsFontType);
            ylabel('feature values',                        ...
                'FontSize', options.labelsFontSize,         ...
                'FontName', options.labelsFontType);

            hold on;
            grid on;
        end
    end
end