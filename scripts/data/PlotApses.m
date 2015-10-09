close all; clear; clc;

load('HaloApses.mat');

figure();
plot3(Apses(1,:), Apses(2,:), Apses(3,:), '^');
hold off; grid on; axis equal;