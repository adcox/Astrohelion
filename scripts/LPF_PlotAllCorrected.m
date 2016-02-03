% Plot all the corrected nodesets
clear; clc; close all;
defineConstants;

dataDir = 'data/LPF_4B_NaturalManifolds/';
files = dir([dataDir, '*_minDV.mat']);
bcSys = bcr4bpR_getSysParam('sun', 'earth', 'moon');
L = bcSys.charL;

fontSize = 18;
lineWidth = 2;

%% XY Plot
figure(); hold on;
for i = 1:(length(files)-1)
    nodes = load([dataDir, files(i).name]);
    orb = bcr4bpR_MS2OrbStruct(nodes);
    plot(orb.State(:,1)*L, orb.State(:,2)*L, 'LineWidth', lineWidth);
    plot(nodes.Nodes(9,1)*L, nodes.Nodes(9,2)*L, 'mv', 'MarkerFaceColor', 'm', 'MarkerSize', 10);
    plot(nodes.Nodes(15,1)*L, nodes.Nodes(15,2)*L, 'mv', 'MarkerFaceColor', 'm', 'MarkerSize', 10);
end
axis equal; limits = axis;
h_env = bcr4bpR_plotEnv('sun', 'earth', 'moon', 'plotsp', 'on', 'scale', L,...
    't0', orb.Time(1), 'tf', orb.Time(end), 'linewidth', lineWidth);
hold off; grid on; axis(limits);
xlabel('x, km');
ylabel('y, km');
set(gca, 'FontSize', fontSize, 'FontWeight', fontWeight);

%% XZ Plot
figure(); hold on;
for i = 1:(length(files)-1)
    nodes = load([dataDir, files(i).name]);
    orb = bcr4bpR_MS2OrbStruct(nodes);
    plot(orb.State(:,1)*L, orb.State(:,3)*L, 'LineWidth', lineWidth);
    plot(nodes.Nodes(9,1)*L, nodes.Nodes(9,3)*L, 'mv', 'MarkerFaceColor', 'm', 'MarkerSize', 10);
    plot(nodes.Nodes(15,1)*L, nodes.Nodes(15,3)*L, 'mv', 'MarkerFaceColor', 'm', 'MarkerSize', 10);
end
axis equal; limits = axis;
plot(h_env.P3Pos(:,1), h_env.P3Pos(:,3), 'Color', col_gray, 'LineWidth', lineWidth);
plot(h_env.SPNoP3(1), h_env.SPNoP3(3), 'r*', 'MarkerSize', 12);
plot(h_env.Handles(1).XData(1), h_env.Handles(1).ZData(1), 'k*');
hold off; grid on; axis(limits);
xlabel('x, km');
ylabel('z, km');
set(gca, 'FontSize', fontSize, 'FontWeight', fontWeight);