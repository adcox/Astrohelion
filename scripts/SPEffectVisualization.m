% Plot SP location perturbations from SE position due to different bodies
clear; clc; close all;
defineConstants;

lineWidth = 2;
fontSize = 20;
markerSize = 8;
load data/SPEffects.mat;
data(1,:) = (data(1,:) - data(1,1))/3600/24; % scale time

numPts = 50;
markedPts = 1:size(data,2)/numPts:size(data,2);
colors = lines(size(data,2));

figure('Position', [0,0,900,600]); hold on;
% figure(); hold on;
plot(data(1,markedPts), data(2,markedPts), '>', 'color', colors(1,:), 'MarkerSize', markerSize); % Moon
plot(data(1,markedPts), data(3,markedPts), '^', 'color', colors(2,:), 'MarkerSize', markerSize); % Mars
plot(data(1,markedPts), data(4,markedPts), 'h', 'color', colors(3,:), 'MarkerSize', markerSize); % Jupiter
plot(data(1,markedPts), data(5,markedPts), 'v', 'color', colors(4,:), 'MarkerSize', markerSize); % Saturn
plot(data(1,markedPts), data(6,markedPts), 'x', 'color', colors(5,:), 'MarkerSize', markerSize); % Venus
plot(data(1,markedPts), data(7,markedPts), 's', 'color', colors(6,:), 'MarkerSize', markerSize); % Mercury
plot(data(1,markedPts), data(8,markedPts), '*', 'color', colors(7,:), 'MarkerSize', markerSize); % Uranus
plot(data(1,markedPts), data(9,markedPts), 'p', 'color', colors(8,:), 'MarkerSize', markerSize); % Neptune

plot(data(1,:), data(2,:), 'color', colors(1,:), 'LineWidth', lineWidth); % Moon
plot(data(1,:), data(3,:), 'color', colors(2,:), 'LineWidth', lineWidth); % Mars
plot(data(1,:), data(4,:), 'color', colors(3,:), 'LineWidth', lineWidth); % Jupiter
plot(data(1,:), data(5,:), 'color', colors(4,:), 'LineWidth', lineWidth); % Saturn
plot(data(1,:), data(6,:), 'color', colors(5,:), 'LineWidth', lineWidth); % Venus
plot(data(1,:), data(7,:), 'color', colors(6,:), 'LineWidth', lineWidth); % Mercury
plot(data(1,:), data(8,:), 'color', colors(7,:), 'LineWidth', lineWidth); % Uranus
plot(data(1,:), data(9,:), 'color', colors(8,:), 'LineWidth', lineWidth); % Neptune
hold off; grid on;
xlabel('Time since 2015/01/01 00:00:00.00 UTC, Days');
ylabel('SP Perturbation Mag. from SE Location, km');
title('SP Location Perturbations Due to Celestial Bodies');
legend('Moon', 'Mars', 'Jupiter', 'Saturn', 'Venus', 'Mercury', 'Uranus', 'Neptune');
set(gca, 'FontSize', fontSize, 'FontWeight', fontWeight);
set(gca, 'yscale', 'log');