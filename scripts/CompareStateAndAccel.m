% Compare detailed state and accel info between Matlab and C++
%
%   To ensure that each state matches exactly, you have to add this code to
%   the bcr4bpR_simNL script after the call to ode113 and before extracting
%   the state and time from data_first:
%
%       o = load('BCTraj.mat');
%       data_first.y = deval(data_first, o.Time);
%       data_first.x = o.Time.';
%
%   Author: Andrew Cox
%   Version: Jan 12, 2015
clear; clc; close all;

bcSys = bcr4bpR_getSysParam('sun', 'earth', 'moon');
o = load('BCTraj.mat');

orbMat = bcr4bpR_simNL(o.State(1,:), 'sun', 'earth', 'moon',...
    't0', o.Time(1), 'tf', o.Time(end) - o.Time(1));

tScalar = bcSys.charT/24/3600;
o.Time = (o.Time - o.Time(1))*tScalar ;
ix = 1:1000;

avgState = mean([abs(orbMat.State), abs(orbMat.Accel)]);

%%
figure(); hold on;
plot(o.Time(ix), (o.State(ix,1) - orbMat.State(ix,1))/avgState(1), '-*');
plot(o.Time(ix), (o.State(ix,2) - orbMat.State(ix,2))/avgState(2), '-*');
plot(o.Time(ix), (o.State(ix,3) - orbMat.State(ix,3))/avgState(3), '-*');
legend('P_x', 'P_y', 'P_z');
hold off; grid on;
xlabel('Epoch, days');
ylabel('Relative Position Difference, non-dim');
title('C++ Pos - Matlab Pos');
set(gca, 'FontSize', 14, 'FontWeight', 'bold');

figure(); hold on;
plot(o.Time(ix), (o.State(ix,4) - orbMat.State(ix,4))/avgState(4), '-*');
plot(o.Time(ix), (o.State(ix,5) - orbMat.State(ix,5))/avgState(5), '-*');
plot(o.Time(ix), (o.State(ix,6) - orbMat.State(ix,6))/avgState(6), '-*');
legend('V_x', 'V_y', 'V_z');
hold off; grid on;
xlabel('Epoch, days');
ylabel('Relative Velocity Difference, non-dim');
title('C++ Vel - Matlab Vel');
set(gca, 'FontSize', 14, 'FontWeight', 'bold');

figure(); hold on;
plot(o.Time(ix), (o.Accel(ix,1) - orbMat.Accel(ix,1))/avgState(7), '-*');
plot(o.Time(ix), (o.Accel(ix,2) - orbMat.Accel(ix,2))/avgState(8), '-*');
plot(o.Time(ix), (o.Accel(ix,3) - orbMat.Accel(ix,3))/avgState(9), '-*');
legend('A_x', 'A_y', 'A_z');
hold off; grid on;
xlabel('Epoch, days');
ylabel('Relative Acceleration Difference, non-dim');
title('C++ Accel - Matlab Accel');
set(gca, 'FontSize', 14, 'FontWeight', 'bold');
