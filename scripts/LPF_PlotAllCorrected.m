% Plot all the corrected nodesets
clear; clc; close all;
defineConstants;

dataDir = 'data/LPF_4B_NaturalManifolds/';
files = dir([dataDir, '*_minDV.mat']);
bcSys = bcr4bpR_getSysParam('sun', 'earth', 'moon');
L = bcSys.charL;

figure(); hold on;
for i = 1:length(files)
    nodes = load([dataDir, files(i).name]);
    orb = bcr4bpR_MS2OrbStruct(nodes);
    plot3(orb.State(:,1)*L, orb.State(:,2)*L, orb.State(:,3)*L, 'LineWidth', lineWidth);
end
axis equal; limits = axis;
bcr4bpR_plotEnv('sun', 'earth', 'moon', 'plotsp', 'on', 'scale', L,...
    't0', orb.Time(1), 'tf', orb.Time(end));
hold off; grid on; axis(limits);
xlabel('x, km');
ylabel('y, km');
zlabel('z, km');
set(gca, 'FontSize', fontSize, 'FontWeight', fontWeight);