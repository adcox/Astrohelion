%% Plot all the arcs from the LPF manifold propagation script
%
%   Promising Candidates
%   9, 11, 12, 57, 58, 67, 78, 84       - Nails a lunar flyby to capture in the EM system
%   51-54       - Conventional transfers, all equally promising
%   59-65, 73, 80-83, 85, 87  - Loops (large) in and out of EM system; may
%       be easier to correct
%   Other arcs
%   21-50   - Nice ribbon of manifolds that flyby and depart

clear; clc; close all;
defineConstants;
colors = parula(10);
fontSize = 18;
lineWidth = 2;

bcSys = bcr4bpR_getSysParam('sun', 'earth', 'moon');
dataDir = 'LPF_4B_NaturalManifolds/';
files = dir([dataDir, 'Traj*_SEM.mat']);
ix0 = 1000; % start plotting data partway through to avoid some of the halo

h_fig = figure(); hold on;
% for i = 1:length(files)
% for i = [9, 11, 12, 57, 58, 67, 78, 84]
% for i = [51:54]
for i = [59:65, 74, 80:83, 85, 87]
% for i = 82;
    o = load([dataDir, files(i).name]);
    % Plot all same color
%     plot3(o.State(ix0:end,1)*bcSys.charL, o.State(ix0:end,2)*bcSys.charL,...
%         o.State(ix0:end,3)*bcSys.charL, 'Color', colors(3,:), 'LineWidth', lineWidth);

    % Plot all different colors
    plot3(o.State(ix0:end,1)*bcSys.charL, o.State(ix0:end,2)*bcSys.charL,...
        o.State(ix0:end,3)*bcSys.charL, 'LineWidth', lineWidth);
    
    axis equal;
end
grid on; limits = axis;
bcr4bpR_plotEnv('Sun', 'Earth', 'Moon', 't0', o.Time(1), 'tf', o.Time(end),...
    'scale', bcSys.charL, 'plotsp', 'on', 'LineWidth', lineWidth);
hold off; axis(limits);
xlabel('x, km');
ylabel('y, km');
zlabel('z, km');
h = findobj(gca, 'type', 'surface');
delete(h);
set(gca, 'FontSize', fontSize, 'FontWeight', fontWeight);



%% Detailed Notes
%
%   9 - Captures VERY well but has significant out-of-plane components. If
%       I can control the flyby, I might be able to tame the out-of-plane bits. 
%   67 - Cool flow pattern that loops around the 3D lunar orbit; probably
%       too out of plane for good SP encounters
%
%   81, 82, 83, 87 - Has very close passes - COULD BE VERY GOOD CANDIDATE
%   84 - Several flybys look close to the right in-plane amplitude