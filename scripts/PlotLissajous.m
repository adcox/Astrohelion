% Plot both the linear and non-linear Lissajous trajectories together
clear; clc; close all;
defineConstants;

liss_lin = load('data/SE_L1_Liss_linear.mat');
liss = load('data/SE_L1_Liss_nonlin.mat');
[mu3B, charT, charL, charM] = cr3bp_getSysParam(liss.P1, liss.P2);

figure(); hold on;
plot3(liss_lin.State(:,1)*charL, liss_lin.State(:,2)*charL,...
    liss_lin.State(:,3)*charL, 'LineWidth', lineWidth);
plot3(liss.State(:,1)*charL, liss.State(:,2)*charL, liss.State(:,3)*charL,...
    'linewidth', lineWidth);
axis equal; limits = axis;
cr3bp_plotEnv(liss.P1, liss.P2);
hold off; axis(limits); grid on;
xlabel('x, km', 'FontSize', fontSize, 'fontWeight', fontWeight);
ylabel('y, km', 'FontSize', fontSize, 'fontWeight', fontWeight);
zlabel('z, km', 'FontSize', fontSize, 'fontWeight', fontWeight);
legend('Linearization', 'Non-Linear');
title('Linear & Non-Linear Lissajous', 'FontSize', fontSize, 'fontWeight', fontWeight);
set(gca, 'FontSize', fontSize, 'FontWeight', fontWeight);