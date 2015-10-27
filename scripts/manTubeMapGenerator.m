% Check the grid of ICs to make sure they fall within the contours
clear; clc; close all;
defineConstants;

load data/manifoldTubeGridding/crossings_C3.150.mat;     % map crossings
load  data/manifoldTubeGridding/gridIC_C3.150.mat;       % Gridded ICs
% load data/apses_C3.020.mat;         % Recorded Apses

tubeMarker = '-k';
gridMarker = 'ob';
gridSize = 3;
gridFill = 'b';

% Add an extra row to gridICs
gridICs = [gridICs; zeros(1, size(gridICs,2))];

count = 1;
for i = 1:size(gridICs,2)
    [I,J] = find(apses(1,count:end) ==  i);
    gridICs(7,i) = length(I);
    count = count + length(I);
end

% Fill that extra row with the number of apses
%% Choosing Z from y
figure(); hold on;
% plot(gridICs(1,:), gridICs(3,:), gridMarker, 'MarkerFaceColor', gridFill, 'MarkerSize', gridSize);
scatter(gridICs(1,:), gridICs(3,:), 3, gridICs(7,:));
plot(crossings(1,:), crossings(3,:), tubeMarker);
hold off; grid on;
xlabel('x, non-dim');
ylabel('z, non-dim');
title('x-z Contour');
set(gca, 'FontSize', fontSize, 'FontWeight', fontWeight);

%% Choosing Vx from Y and Z
figure('Position', [50,50,1200,500]);
subplot(1,2,1); hold on;
plot(gridICs(1,:), gridICs(4,:), gridMarker, 'MarkerFaceColor', gridFill, 'MarkerSize', gridSize);
plot(crossings(1,:), crossings(4,:), tubeMarker);
hold off; grid on;
xlabel('x, non-dim');
ylabel('v_x, non-dim');
title('x-v_x Contour');
set(gca, 'FontSize', fontSize, 'FontWeight', fontWeight);

subplot(1,2,2); hold on;
plot(gridICs(3,:), gridICs(4,:), gridMarker, 'MarkerFaceColor', gridFill, 'MarkerSize', gridSize);
plot(crossings(3,:), crossings(4,:), tubeMarker);
hold off; grid on;
xlabel('z, non-dim');
ylabel('v_x, non-dim');
title('z-v_x Contour');
set(gca, 'FontSize', fontSize, 'FontWeight', fontWeight);

%% Choosing Vy from Y, Z, and Vx
figure('Position', [50,50,1200,900]);
subplot(2,2,1); hold on;
plot(gridICs(1,:), gridICs(5,:), gridMarker, 'MarkerFaceColor', gridFill, 'MarkerSize', gridSize);
plot(crossings(1,:), crossings(5,:), tubeMarker);
hold off; grid on;
xlabel('x, non-dim');
ylabel('v_y, non-dim');
title('x-v_y Contour');
set(gca, 'FontSize', fontSize, 'FontWeight', fontWeight);

subplot(2,2,2); hold on;
plot(gridICs(3,:), gridICs(5,:), gridMarker, 'MarkerFaceColor', gridFill, 'MarkerSize', gridSize);
plot(crossings(3,:), crossings(5,:), tubeMarker);
hold off; grid on;
xlabel('z, non-dim');
ylabel('v_y, non-dim');
title('z-v_y Contour');
set(gca, 'FontSize', fontSize, 'FontWeight', fontWeight);

subplot(2,2,3); hold on;
plot(gridICs(4,:), gridICs(5,:), gridMarker, 'MarkerFaceColor', gridFill, 'MarkerSize', gridSize);
plot(crossings(4,:), crossings(5,:), tubeMarker);
hold off; grid on;
xlabel('v_x, non-dim');
ylabel('v_y, non-dim');
title('v_x-v_y Contour');
set(gca, 'FontSize', fontSize, 'FontWeight', fontWeight);

%% Choosing Vz from C and z
figure('Position', [50,50,1200,900]);
subplot(2,2,1); hold on;
plot(gridICs(1,:), gridICs(6,:), gridMarker, 'MarkerFaceColor', gridFill, 'MarkerSize', gridSize);
plot(crossings(1,:), crossings(6,:), tubeMarker);
hold off; grid on;
xlabel('x, non-dim');
ylabel('v_z, non-dim');
title('x-v_z Contour');
set(gca, 'FontSize', fontSize, 'FontWeight', fontWeight);

subplot(2,2,2); hold on;
plot(gridICs(3,:), gridICs(6,:), gridMarker, 'MarkerFaceColor', gridFill, 'MarkerSize', gridSize);
plot(crossings(3,:), crossings(6,:), tubeMarker);
hold off; grid on;
xlabel('z, non-dim');
ylabel('v_z, non-dim');
title('z-v_z Contour');
set(gca, 'FontSize', fontSize, 'FontWeight', fontWeight);

subplot(2,2,3); hold on;
plot(gridICs(4,:), gridICs(6,:), gridMarker, 'MarkerFaceColor', gridFill, 'MarkerSize', gridSize);
plot(crossings(4,:), crossings(6,:), tubeMarker);
hold off; grid on;
xlabel('v_x, non-dim');
ylabel('v_z, non-dim');
title('v_x-v_z Contour');
set(gca, 'FontSize', fontSize, 'FontWeight', fontWeight);

subplot(2,2,4); hold on;
plot(gridICs(5,:), gridICs(6,:), gridMarker, 'MarkerFaceColor', gridFill, 'MarkerSize', gridSize);
plot(crossings(5,:), crossings(6,:), tubeMarker);
hold off; grid on;
xlabel('v_y, non-dim');
ylabel('v_z, non-dim');
title('v_y-v_z Contour');
set(gca, 'FontSize', fontSize, 'FontWeight', fontWeight);
