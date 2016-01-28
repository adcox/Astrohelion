% Plot a loopy trajectory
clear; close all; clc;
defineConstants;
colors = parula(10);
lineWidth = 2;
fontSize = 18;

% traj = load('data/Traj066.mat');
traj = load('data/Traj008.mat');
ix0 = 1000; % don't plot the first points to skip the corrected halo

bcSys = bcr4bpR_getSysParam('sun', 'earth', 'moon');
L = bcSys.charL;

figure(1); hold on;
plot(traj.State(ix0:end,1)*L, traj.State(ix0:end,2)*L,...
    'Color', colors(3,:), 'LineWidth', lineWidth);
axis equal; limits = axis;
envData = bcr4bpR_plotEnv('sun', 'earth', 'moon', 't0', traj.Time(1),...
    'tf', traj.Time(end), 'theta0', traj.Theta0, 'phi0', traj.Phi0,...
    'gamma', traj.Gamma, 'plotsp', 'on', 'scale', L, 'LineWidth', lineWidth);
hold off; grid on;
axis(limits + 0.1*[-diff(limits(1:2)), diff(limits(1:2)), -diff(limits(3:4)), diff(limits(3:4))]);
xlabel('x, km');
ylabel('y, km');
set(gca, 'FontSize', fontSize, 'FontWeight', fontWeight);
h = findobj(gca, 'Type', 'Surface');
delete(h);

%%
figure(2); hold on;
plot(envData.P3Pos(:,1), envData.P3Pos(:,3), 'Color', col_gray, 'LineWidth', lineWidth);
plot(0, 0, 'ok', 'MarkerFaceColor', 'k');
plot(envData.SPPos(end,1), envData.SPPos(end,3), 'r*', 'MarkerSize', 15);
plot(traj.State(ix0:end,1)*L, traj.State(ix0:end,3)*L, 'Color', colors(3,:), 'LineWidth', lineWidth);
axis equal; hold off; grid on;
xlabel('x, km');
ylabel('z, km');
set(gca, 'FontSize', fontSize, 'FontWeight', fontWeight);