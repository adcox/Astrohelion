function [] = plotMap(datafile, varargin)
% PLOTMAP Plot a Poincare Map
%
%   plotMap(datafile) plots a map of the data given in the specified
%       file.
%
%   plotMap(datafile, ...) does the same, but with additional parameters
%
%   PARAMETERS
%
%       'Highlight'     -   The orbit number you'd like to highlight
%
%       'Labeled'       -   The number of returns on the highlighed orbit
%                           you would like to label
%
%       'Range'         -   The axes range you'd like to view. Should be a 
%                           proper input to the AXIS command. Default value 
%                           is [-1 1 -5 5]
%
%       'SaveStr'       -   Name of the file that will be saved. Default is
%                           'map'.
%
%   See also AXIS

%% Default options & Input Handling
def_range = [-2 1.5 -5 5];
def_saveStr = 'map';
def_debug = 'Off';
def_highlight = 0;
def_labelReturns = 0;

isFourVector = @(x) length(x) == 4 && isnumeric(x);
isPosNum = @(x) isscalar(x) && x > 0;

p = inputParser;

addRequired(p,'datafile', @isstr);

addParameter(p, 'Range', def_range, isFourVector);
addParameter(p, 'saveStr', def_saveStr, @isstr);
addParameter(p, 'ShowDebug', def_debug, @isstr);
addParameter(p, 'Highlight', def_highlight, isPosNum);
addParameter(p, 'Labeled', def_labelReturns, isPosNum);

parse(p, datafile, varargin{:});

range = p.Results.Range;
saveStr = p.Results.saveStr;
showDebug = strcmpi(p.Results.ShowDebug, 'On');
highlight = p.Results.Highlight;
labelReturns = p.Results.Labeled;

load(datafile);
defineConstants;

% Pick a dot size based on range
defMaxDim = max([abs(def_range(2)-def_range(1)), abs(def_range(4)-def_range(3))]);
maxDim = max([abs(range(2)-range(1)), abs(range(4)-range(3))]);
k = defMaxDim/maxDim;
if(k <= 1)
    dotSize = 2.5;
else
    dotSize = 2.5 + 0.5*(k-1);
end

if(highlight > 0)
    orbsToHL = find(orbitNums == highlight);
end

%% Plotting
f_x_xdot_map= figure(); hold on;
H_points = plot(allReturns(:,1), allReturns(:,4), '.', 'MarkerSize', dotSize);
H_bound = plot(constCPts(:,1), constCPts(:,2), 'Color', col_gray, 'LineWidth', lineWidth);
plot(constCPts(:,1), -constCPts(:,2), 'Color', col_gray, 'LineWidth', lineWidth); 

if(highlight > 0)
    plot(allReturns(orbsToHL,1), allReturns(orbsToHL,4), 'm.', 'MarkerSize', 1.5*dotSize);
    
    if(labelReturns > 0)
        plot(allReturns(orbsToHL(1:labelReturns),1), allReturns(orbsToHL(1:labelReturns),4),...
            'k.', 'MarkerSize', 2*dotSize);
        for n = 1:labelReturns
            text(allReturns(orbsToHL(n),1), allReturns(orbsToHL(n),4),...
                sprintf('%d', n),...
                'horizontalAlignment', 'center',...
                'verticalAlignment', 'bottom', ...
                'FontSize', 6/7*fontSize, 'FontWeight', fontWeight);
        end
    end
end

hold off; grid on;
set(gca, 'YLim', range(3:4), 'XLim', range(1:2));
xlabel('x, non-dim');
ylabel('x-dot, non-dim');
title(sprintf('Poincare Map for C = %.4f', C));
set(gca, 'FontSize', 18, 'FontWeight', 'bold');
% saveas(f_x_xdot_map, sprintf('plots/%s', saveStr), 'png');
    
end