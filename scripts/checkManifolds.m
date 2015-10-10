% Plot manifolds
clear; clc; close all;
defineConstants;

% load data/ManCross_IC.mat;
load data/crossings_C3.030.mat;

[mu3B, charT,charL,charM] = cr3bp_getSysParam('earth', 'moon');
LPts = cr3bp_getEquilibPts(mu3B);

f_3d = figure(); hold on;

% draw line at proposed hyperplane
y = LPts.L2(2) + 0.3;
planePoints = [2 y 0.5; 1.8 y -0.5; 0.8 y -0.5; 0.8 y 0.5];
fill3(planePoints(:,1), planePoints(:,2), planePoints(:,3), 'g');
alpha(0.3);

for i = 1:10:max(size(crossings))
    fprintf('Simulating %03d/%03d...\n', i, max(size(crossings)));
    o = cr3bp_simNL(crossings(1:6,i), 'earth', 'moon', 'tf', 2*pi);
    
    plot3(o.State(:,1), o.State(:,2), o.State(:,3), 'color', col_green, ...
        'linewidth', lineWidth);
end
view(3); axis equal;
plot3(crossings(1,:), crossings(2,:), crossings(3,:), '*k');


limits = [0.8, 1.8, -0.5, 0.5, -0.5, 0.5];
cr3bp_plotEnv('earth', 'moon');
hold off; grid on; axis(limits);
xlabel('x, non-dim', 'FontSize', fontSize, 'fontweight', fontWeight);
ylabel('y, non-dim',  'FontSize', fontSize, 'fontweight', fontWeight);
zlabel('z, non-dim',  'FontSize', fontSize, 'fontweight', fontWeight);
title('L2 S-Halo Stable (+) Manifolds (Rev Time)',  'FontSize', fontSize,...
    'fontweight', fontWeight);
set(gca, 'FontSize', fontSize, 'fontweight', fontWeight);

%% Contours
figure();
plot(crossings(1,:), crossings(3,:), '-*');
axis equal;
xlabel('x');
ylabel('z');

figure();
plot(crossings(1,:), crossings(5,:), '-*');
xlabel('x');
ylabel('v_y');

figure();
plot(crossings(1,:), crossings(6,:), '-*');
xlabel('x');
ylabel('v_z');