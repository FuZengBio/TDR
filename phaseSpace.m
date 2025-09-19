% Code orignially from Mark Churchland, adapted by Adam Zaidel
% For making publication quality rosette plots
% useage:
%       phaseSpace(Projection, Summary)
%       phaseSpace(Projection, Summary, params)
%
% You can limit the conditions you plot either by
%   1) including only some entries of Projection (e.g., Projection(1:2))
%   or 2) by using params.conds2plot to restrict (e.g., params.conds2plot = 1:2);
%   In the former case, scaling etc is based just on the passed points.
%   In the latter scaling is based on all the points.
%
% outputs are [colorStruct, hf, haxP, vaxP] = phaseSpace(Projection, Summary)
%
%   'colorStruct' is a structure (one per dataset) of cells of linecolors (one per condition)
%   that you might wish to pass to another function (e.g., one that plots the rosette PSTH's
%   or the hand trajectories).
%
%    'hf' is the fig#, 'haxP' and 'vaxP' are the axis parameters.
%
%   *** ALL of the outputs, except the last, pertain only the LAST graph plotted. ***
%
% rosetteData comes from multiRosetteScript
%
% params can have the following fields:
%       .times    These override the default times (those corresponding to scores; e.g., the orginal
%               times that were used to build the space.  If empty, the defaults are used.
%                 Note that zero is movement onset.  Thus .times will probably start strongly
%               negatively.
%                 Thus, you might pass it -1550:10:150 if you wanted to start way back at the
%               beginning.
%                 A nice feature is that only those times that match times in 'scoresExtraTime' are
%               used.  Thus, if you pass -100000:1000000 things will still work fine.
%               Scalings are based on all times in Projection.projAllTimes.  This is nice for movies
%               as the scaling won't change as a function of the times you plot.
%
%       .planes2plot   list of the jPC planes you want plotted.  Default is [1].  [1,2] would also be reasonable.
%
%       .arrowSize        The default is 5
%       .arrowGain        FOR MOVIES: sets velocity dependence of arrow size (0 to not grow when faster).
%       .plotPlanEllipse  Controls whether the ellipse is plotted.  The default is 'true'
%       .useAxes          whether axes should be plotted.  Default is 'true'
%       .useLabel         whether to label with dataset. Default is 'true'
%       .planMarkerSize   size of the plan dot.  Default is 6.
%       .lineWidth        width of the trajectories.  Default is 0.85.
%       .arrowMinVel      minimum velocity for plotting an arrow.
%       .rankType         default is 'eig', but you can override with 'varCapt'.  The first plane will
%                         then be the jPC plane that captured the most variance (often associated with the
%                         largest eigenvalue but not always
%       .conds2plot       which conditions to plot (scalings will still be based on all the conds in 'Projection')
%       .substPCs         0 = use jPCA projections; 1 = use rawPC projections; 2=use Targeted PCS!
%       .crossCondMean    if present and == 1, plot the cross condition mean in cyan.
%       .reusePlot        if present and == 1, do cla then reuse the plot
%       .dataRanges       normally this is set automatically, but you can decide yourself what the
%                         range should be.  You should supply one entry per plane to be plotted.  You can also supply
%                         just the first, and then the defaults will be used after that.
%
function [colorStruct, haxP, vaxP, ph] = phaseSpace(Projection, Summary, params)


%% some basic parameters
axLimScale = 1.1;
axisSeparation = 0.10;  % separated by 20% of the maximum excursion (may need more if plotting future times, which aren't used to compute farthestLeft or farthestDown)

numPlanes = length(Summary.varCaptEachPlanePC);  % total number of planes provided (may only plot a subset)

%% set defaults and override if 'params' is included as an argument
plot_steps=0;
% allows for the use of times other than the original ones that correspond to 'scores' (those that
% were used to create the projection and do the analysis)
overrideTimes = [];
if exist('params', 'var') && isfield(params,'times')
    overrideTimes = params.times;
end

if exist('params', 'var') && isfield(params,'RegNames')
    RegNames=nameREG(params.RegNames);
end

arrowSize = 5;
if exist('params', 'var') && isfield(params,'arrowSize')
    arrowSize = params.arrowSize;
end

arrowGain = 0;
if exist('params', 'var') && isfield(params,'arrowGain')
    arrowGain = params.arrowGain;
end

% Default is we plot the ellipse if we have 6 or more conditions
% NOTE: we still plot it even if asked to only plot one cond, so long as we HAVE more than 6 to
% build the ellipse off.  It matters whether length(Projection) >= 6, not whether conds2plot >= 6
if length(Projection) >= 6
    plotPlanEllipse = true;
else
    plotPlanEllipse = false;
    axLimScale = 1.3*axLimScale;
end
if exist('params', 'var') && isfield(params,'plotPlanEllipse')
    plotPlanEllipse = params.plotPlanEllipse;
end

useAxes = true;
if exist('params', 'var') && isfield(params,'useAxes')
    useAxes = params.useAxes;
end

ArrowTime = false;
if exist('params', 'var') && isfield(params,'ArrowTime')
    ArrowTime = params.ArrowTime;
end

useLabel = true;
if exist('params', 'var') && isfield(params,'useLabel')
    useLabel = params.useLabel;
end

planMarkerSize = 6;
if exist('params', 'var') && isfield(params,'planMarkerSize')
    planMarkerSize = params.planMarkerSize;
end

lineWidth = 0.85;
if exist('params', 'var') && isfield(params,'lineWidth')
    lineWidth = params.lineWidth;
end

arrowMinVel = [];
if exist('params', 'var') && isfield(params,'arrowMinVel')
    arrowMinVel = params.arrowMinVel;
end

planes2plot = 1;  % this is a list of which planes to plot
if exist('params', 'var') && isfield(params,'planes2plot')
    planes2plot = params.planes2plot;
end

%default?
if exist('params', 'var') && isfield(params,'params2plot')
    params2plot = params.params2plot;
end

rankType = 'eig';
if exist('params', 'var') && isfield(params,'rankType')
    rankType = params.rankType;
end

reusePlot = 0;
if exist('params', 'var') && isfield(params,'reusePlot')
    reusePlot = params.reusePlot;
end
if length(planes2plot) > 1, reusePlot = 0; end  % cant reuse if we are plotting more than one thing.

numConds = length(Projection);
conds2plot = 1:numConds;
if exist('params', 'var') && isfield(params,'conds2plot')
    if ~strcmp(params.conds2plot,'all')
        conds2plot = params.conds2plot;
    end
end

% if asked, plot the cross-condition mean on the same plot.
crossCondMean = false;
if exist('params', 'var') && isfield(params,'crossCondMean')
    crossCondMean = params.crossCondMean;
end

% If asked, substitue the raw PC projections.
substPCs = 0;
if exist('params', 'var') && isfield(params,'substPCs') && params.substPCs
    substPCs = params.substPCs;
    if substPCs==1
        % just overwrite
        for c = 1:numConds
            Projection(c).proj = Projection(c).tradPCAproj;
            Projection(c).projAllTimes = Projection(c).tradPCAprojAllTimes;
        end
        Summary.varCaptEachPlane = sum(reshape(Summary.varCaptEachPC,2,numPlanes));
    elseif substPCs==2
        for c = 1:numConds
            Projection(c).proj = Projection(c).TargetedPCAproj;
            Projection(c).projAllTimes = Projection(c).TargetedPCAprojAllTimes;
        end
        %TBD currently I only use 2 params (params2plot)
        Summary.varCaptEachPlane = [sum(Summary.varCaptEachTargetedPC(params2plot))];
    end
end

if strcmp(rankType, 'varCapt')  && substPCs == 0  % we WONT reorder if they were PCs
    [~, sortIndices] = sort(Summary.varCaptEachPlane,'descend');
    planes2plot_Orig = planes2plot;  % keep this so we can label the plane appropriately
    planes2plot = sortIndices(planes2plot);  % get the asked for planes, but by var accounted for rather than eigenvalue
end

% the range of the data will set the size of the plot unless you manually override
dataRanges = max(abs(vertcat(Projection.proj)));
if substPCs==2
    %TBD currently I only use 2 params (params2plot)
    dataRanges = max(dataRanges(params2plot));
else
    dataRanges = max(reshape(dataRanges,2,numPlanes));  % one range per plane
end

if exist('params', 'var') && isfield(params,'dataRanges')
    for i = 1:length(params.dataRanges)
        dataRanges(i) = params.dataRanges(i);  % only override those values that are specified
    end
end



arrowEdgeColor = 'k';

for pindex = 1:length(planes2plot)
    
    % get some useful indices
    plane = planes2plot(pindex);  % which plane to plot
    
    if exist('params2plot', 'var') %added by AZ
        d1=params2plot(1);
        d2=params2plot(2);
    else %the original way of choosing dimensions
        d2 = 2*plane;  % indices into the dimensions
        d1 = d2-1;
    end
    
    % set the limits of the figure
    axLim = axLimScale * dataRanges(plane) * [-1 1 -1 1];
    if axLimScale * dataRanges(plane)>50 %if not normalized
        axisLength=50;
    elseif axLimScale * dataRanges(plane)>25 %if not normalized
        axisLength=25;
    elseif axLimScale * dataRanges(plane)>10 %if not normalized
        axisLength=10;
    else %normalized
        axisLength = 0.5;
    end
    
    
    for c = 1:numConds
        % Always taken from the 1st element (NOT the first that will be plotted)
        % This way the ellipse doesn't depend on which times you choose to plot
        planData(c,:) = Projection(c).proj(1,[d1,d2]);
    end
    % need this for ellipse plotting
    if numConds > 1   % used to be only if plotPlanEllipse==1.  Now we always use the same scaling regardless of whether we use the plan ellipse.
        ellipseRadii = 2*var(planData).^0.5;  % we may plot an ellipse for the plan activity
        
        % these will be altered further below based on how far the data extends left and down
        farthestLeft = -ellipseRadii(1);  % used figure out how far the axes need to be offset
        farthestRight = ellipseRadii(1);  % used figure out how far the axes need to be offset
        farthestDown = -ellipseRadii(2);  % used figure out how far the axes need to be offset
        farthestUp = ellipseRadii(2);  % used figure out how far the axes need to be offset
    else
        temp = vertcat(Projection.proj);
        temp = temp(:,[d1,d2]);
        farthestLeft = -1.2*max(temp(:));  % used figure out how far the axes need to be offset
        farthestRight = 1.2*max(temp(:));  % used figure out how far the axes need to be offset
        farthestDown = -1.2*max(temp(:));  % used figure out how far the axes need to be offset
        farthestUp = 1.2*max(temp(:));  % used figure out how far the axes need to be offset
    end
    
    %% deal with the color scheme
    
    % ** colors graded based on PLAN STATE
    % These do NOT depend on which times you choose to plot (only on which time is first in Projection.proj).
    htmp = redgreencmap(numConds, 'interpolation', 'linear');
    [~,newColorIndices] = sort(planData(:,1));
    
    htmp(newColorIndices,:) = htmp;
    
    for c = 1:numConds  % cycle through conditions, and assign that condition's color
        lineColor{c} = htmp(c,:);
        arrowFaceColor{c} = htmp(c,:);
        planMarkerColor{c} = htmp(c,:);
        lineStyle{c} = '-'; %default
    end
    
    % override colors if asked
    if exist('params', 'var') && isfield(params,'colors')
        lineColor = params.colors;
        arrowFaceColor = params.colors;
        planMarkerColor = params.colors;
    end
    
    colorStruct(pindex).colors = lineColor;
    
    % override linestyle if asked
    if exist('params', 'var') && isfield(params,'style')
        lineStyle = params.style;
    end
    
    %% Plot the rosette itself
    if reusePlot == 0, blankFigure(axLim); else cla; end
    daspect([1 1 1])
    
    % first deal with the ellipse for the plan variance (we want this under the rest of the data)
    if plotPlanEllipse, circle([0 0], ellipseRadii, 0.6*[1 1 1], 1); end
    
    % cycle through conditions
    for c = 1:numConds
        
        if isempty(overrideTimes)  % if we are going with the original times (those that were used to create the projection and do the analysis)
            Times=Projection(c).allTimes;
            P1 = Projection(c).proj(:,d1);
            P2 = Projection(c).proj(:,d2);
        else
            useTimes = ismember(Projection(c).allTimes, overrideTimes);
            Times=overrideTimes;
            P1 = Projection(c).projAllTimes(useTimes,d1);
            P2 = Projection(c).projAllTimes(useTimes,d2);
        end
        
        if ismember(c,conds2plot)
            
            if plot_steps
                ph(c)=plot(P1, P2, 'ko', 'markerSize', planMarkerSize, 'markerFaceColor', planMarkerColor{c}, 'color', lineColor{c}, 'lineWidth', lineWidth, 'linestyle',lineStyle{c});
            else
                ph(c)=plot(P1, P2, 'color', lineColor{c}, 'lineWidth', lineWidth, 'linestyle',lineStyle{c});
                if planMarkerSize>0
                    plot(P1(1), P2(1), 'ko', 'markerSize', planMarkerSize, 'markerFaceColor', planMarkerColor{c});
                end
            end
            
            
            % for arrow, figure out last two points, and (if asked) supress the arrow if velocity is
            % below a threshold.
            penultimatePoint = [P1(end-1), P2(end-1)];
            lastPoint = [P1(end), P2(end)];
            vel = norm(lastPoint - penultimatePoint);
            if isempty(arrowMinVel) || vel > arrowMinVel
                aSize = arrowSize + arrowGain * vel;  % if asked (e.g. for movies) arrow size may grow with vel
                arrowMMC(penultimatePoint, lastPoint, [], aSize, axLim, arrowFaceColor{c}, arrowEdgeColor);
            else
                plot(lastPoint(1), lastPoint(2), 'ko', 'markerSi', arrowSize, 'markerFac', arrowFaceColor{c}, 'markerEdge', arrowEdgeColor);
            end
        end
        
        % axis locations will be based on the original set of times used to make the scores
        % and not on the actual times used.  Here we get the leftmost and bottommost point
        if isfield(Projection, 'projAllTimes')
            farthestLeft = min(farthestLeft, min(Projection(c).projAllTimes(:,d1)));
            farthestRight = max(farthestRight, max(Projection(c).projAllTimes(:,d1)));
            farthestDown = min(farthestDown, min(Projection(c).projAllTimes(:,d2)));
            farthestUp = max(farthestUp, max(Projection(c).projAllTimes(:,d2)));
        else
            farthestLeft = min(farthestLeft, min(Projection(c).proj(:,d1)));
            farthestRight = max(farthestRight, min(Projection(c).proj(:,d1)));
            farthestDown = min(farthestDown, min(Projection(c).proj(:,d2)));
            farthestUp = min(farthestUp, min(Projection(c).proj(:,d2)));
        end
    end
    
    plot(0,0,'b+', 'markerSi', 7.5);  % plot a central cross
    %to see PC vectors
    %figure; plot(Summary.PCs(:,1:2))
    
    
    %% if asked we will also plot the cross condition mean
    if crossCondMean && length(Summary.crossCondMean) > 1
        meanColor = [0 1 1];
        
        if isempty(overrideTimes)  % if we are going with the original times (those that were used to create the projection and do the analysis)
            P1 = Summary.crossCondMean(:,d1);
            P2 = Summary.crossCondMean(:,d2);
        else
            useTimes = ismember(Projection(c).allTimes, overrideTimes);
            P1 = Summary.crossCondMeanAllTimes(useTimes,d1);
            P2 = Summary.crossCondMeanAllTimes(useTimes,d2);
        end
        
        plot(P1, P2, 'color', meanColor, 'lineWidth', 1.2*lineWidth);  % make slightly thicker than for rest of data.
        if planMarkerSize>0
            plot(P1(1), P2(1), 'ko', 'markerSize', planMarkerSize, 'markerFaceColor', meanColor);
        end
        
        % for arrow, figure out last two points, and (if asked) supress the arrow if velocity is
        % below a threshold.
        penultimatePoint = [P1(end-1), P2(end-1)];
        lastPoint = [P1(end), P2(end)];
        vel = norm(lastPoint - penultimatePoint);
        
        aSize = arrowSize + arrowGain * vel;  % if asked (e.g. for movies) arrow size may grow with vel
        arrowMMC(penultimatePoint, lastPoint, [], aSize, axLim, meanColor, arrowEdgeColor);
        
    end
    
    %% make axes
    if useAxes
        clear axisParams;
        
        extraSeparation = axisSeparation*(min(farthestDown,farthestLeft));
        
        % general axis parameters
        axisParams.tickLocations = [-axisLength, 0, axisLength];
        axisParams.longTicks = 0;
        axisParams.fontSize = 12;
        
        % horizontal axis
        axisParams.axisOffset = min(farthestDown + extraSeparation,-axisLength*1.2);
        if substPCs==1
            axisParams.axisLabel = 'Projection onto PC_1 [a.u.]';
        elseif substPCs==2
            axisParams.axisLabel = sprintf('%s (%d%% var)', RegNames{d1},round(100*Summary.varCaptEachTargetedPC(d1)));
            %axisParams.axisLabel = sprintf('%s [a.u.]', RegNames{d1});
        else
            axisParams.axisLabel = 'Projection onto jPC_1 [a.u.]';
        end
        axisParams.axisOrientation = 'h';
        axisParams.axisLabelOffset = axisLength*0.3;
        haxP = AxisMMC(-axisLength, axisLength, axisParams);
        
        % vertical axis
        axisParams.axisOffset = min(farthestLeft + extraSeparation,-axisLength*1.2);
        if substPCs==1
            axisParams.axisLabel = 'Projection onto PC_2 [a.u.]';
        elseif substPCs==2
            axisParams.axisLabel = sprintf('%s (%d%% var)', RegNames{d2},round(100*Summary.varCaptEachTargetedPC(d2)));
            %axisParams.axisLabel = sprintf('%s [a.u.]', RegNames{d2});
        else
            axisParams.axisLabel = 'Projection onto jPC_2 [a.u.]';
        end
        axisParams.axisOrientation = 'v';
        axisParams.axisLabelOffset = axisLength*0.5;
        vaxP = AxisMMC(-axisLength, axisLength, axisParams);        
    end
    %plot arrow time
    if ArrowTime
        text(vaxP.labelx,haxP.labely,sprintf('t=%.1fs',Times(end)/1000),'fontsize',12)
    end
    
    % plot a label at the top
    if substPCs == 0, planeType = 'jPCA'; else planeType = 'PCA'; end
    if useLabel
        if substPCs == 1
            titleText = sprintf('raw PCA plane %d', plane);
        elseif substPCs == 2
            titleText = sprintf('Targeted PCA');
        elseif strcmp(rankType, 'varCapt')
            letters = 'ABCDEFGHIJKLMNOPQRSTUVWXYZ';
            titleText = sprintf('jPCA plane %s (var capt ranked)', letters(planes2plot_Orig(pindex)));
        else
            titleText = sprintf('jPCA plane %d (eigval ranked)', plane);
        end
        %         titleText2 = sprintf('%d%% of var captured', round(100*Summary.varCaptEachPlane(plane)));
        %         title(sprintf('%s\n%s',titleText,titleText2))
        title(sprintf('%s',titleText))
    end
end  % done looping through planes

end  % end of the main function




