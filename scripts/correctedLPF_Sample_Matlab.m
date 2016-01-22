% Attempt to correct the LPF Trajectory using Matlab
clear; clc; close all;
defineConstants;

% Load the data and make a plot
n = load('WorkingNodes.mat');
o = bcr4bpR_MS2OrbStruct(n);

sp1Node = 10;
sp2Node = 15;

% Allow delta V's at beginning of manifold and once on the added arc
velCon = struct('type', 'AllBut', 'data', [sp1Node-4, sp2Node-2]);

bcr4bpR_plotOrbit(o); hold on;
plot3(n.Nodes(:,1), n.Nodes(:,2), n.Nodes(:,3), 'bd');
plot3(n.Nodes(sp1Node,1), n.Nodes(sp1Node,2), n.Nodes(sp1Node,3), 'mo');
plot3(n.Nodes(sp2Node,1), n.Nodes(sp2Node,2), n.Nodes(sp2Node,3), 'mo');
plot3(n.Nodes(velCon.data,1), n.Nodes(velCon.data,2), n.Nodes(velCon.data,3), 'r*');
bcr4bpR_plotEnv('sun', 'earth', 'moon', 't0', n.Epochs(1), 'tf', n.Epochs(end)+n.TOFs(end), 'plotsp', 'on');
f_plot = gcf;
hold off;

%% Correct
% SP Targeting constraints
sp1 = struct('node', sp1Node, 'type', 'sp', 'data', NaN);
sp2 = struct('node', sp2Node, 'type', 'sp', 'data', NaN);

% Fix the first node in space and time to ensure continuity with fixed Halo
fixHalo = struct('type', 'state', 'node', 1, 'data', [n.Nodes(1,:), n.Epochs(1)]);

firstPassData = bcr4bpR_multShoot(n.Nodes, n.TOFs, 'Sun', 'Earth', 'Moon',...
    'Epochs', n.Epochs, 'velcon', velCon, 'tol', 1e-10,...
    'constraints', [sp1, sp2, fixHalo], 'showIts', 'off');
firstPassOrb = bcr4bpR_MS2OrbStruct(firstPassData);

figure(f_plot); hold on;
plot3(firstPassOrb.State(:,1), firstPassOrb.State(:,2), firstPassOrb.State(:,3),...
    'g', 'LineWidth', lineWidth);
plot3(firstPassData.Nodes(:,1), firstPassData.Nodes(:,2), firstPassData.Nodes(:,3),...
    'gd', 'LineWidth', lineWidth);
% bcr4bpR_plotEnv('sun', 'earth', 'moon', firstPassData.Epochs(1), firstPassData.Epochs(2), 'plotsp', 'on');
hold off;