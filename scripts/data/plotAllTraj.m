%% Plot all the arcs from the LPF manifold propagation script
%
%   Promising Candidates
%   9, 11, 12, 57, 58, 67, 78, 84       - Nails a lunar flyby to capture in the EM system
%   51-54       - Conventional transfers
%   59-65, 73, 80-83, 85, 87  - Loops (large) in and out of EM system; large-ish z-comp.
%
%   Other arcs
%   21-50   - Nice ribbon of manifolds that flyby and depart

clear; clc; close all;
defineConstants;
bcSys = bcr4bpR_getSysParam('sun', 'earth', 'moon');

files = dir('LPF_4B_NaturalManifolds/Traj*.mat');

h_fig = figure(); hold on;
% for i = 1:length(files)
for i = 91:100
    o = load(['LPF_4B_NaturalManifolds/', files(i).name]);
    plot3(o.State(:,1)*bcSys.charL, o.State(:,2)*bcSys.charL,...
        o.State(:,3)*bcSys.charL, 'LineWidth', lineWidth);
    axis equal;
end
grid on; limits = axis;
bcr4bpR_plotEnv('Sun', 'Earth', 'Moon', 't0', o.Time(1), 'tf', o.Time(end),...
    'scale', bcSys.charL, 'plotsp', 'on');
hold off; axis(limits);
xlabel('x');
ylabel('y');
zlabel('z');

%% integrate further
o = bcr4bpR_simNL(o.State(1,:), 'sun', 'earth', 'moon',...
    't0', o.Time(1), 'tf', o.Time(end)-o.Time(1));