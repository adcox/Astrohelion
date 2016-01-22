% Plot SP location as function of time, compare to 2nd order polyfit from
% C++
clear; clc; close all;

bcSys = bcr4bpR_getSysParam('sun', 'earth', 'moon');
% T0 = 25.660;
T0 = 10232.85;
tSpan = 1*24*3600/bcSys.charT;
numPts = 100;

load PolyFitCoeff.csv;

T = linspace(T0 - tSpan, T0 + tSpan, numPts);
allSP_pos = zeros(numPts, 3);
approxPos = zeros(numPts, 3);
for i = 1:length(T)
    if(i == 1)
        p = 0;
    else
        p = allSP_pos(i-1,:);
    end
    
    allSP_pos(i,:) = bcr4bpR_locateSP(T(i), p, 0, 0, 5.14*pi/180, bcSys);
    approxPos(i,:) = [T(i)*T(i), T(i), 1]*PolyFitCoeff;
end

%% Attempt modified scaling
t = linspace(-tSpan, tSpan, numPts);
V = [t.'.*t.', t.', ones(numPts, 1)];
C = inv(V.'*V)*V.'*allSP_pos;

linTerms = C(2,:) - 2*C(1,:)*T0;
shiftTerms = C(3,:) + C(1,:)*T0*T0 - C(2,:)*T0;
C0 = C;
C = [C(1,:); linTerms; shiftTerms];

approxPos2 = [T.'.*T.', T.', ones(numPts, 1)]*C;

%% Make plots to check
figure();
subplot(3,1,1); hold on;
plot(T, allSP_pos(:,1), 'LineWidth', 2);
plot(T, approxPos(:,1), '--', 'LineWidth', 2);
% plot(T, approxPos2(:,1), '-.', 'LineWidth', 2);
hold off; grid on;
ylabel('sp_x, non-dim');
% legend('True', 'C++ Approx', 'Centered Approx');
legend('True', 'Centered Approx');
title('Saddle Point Position over Time');
set(gca, 'FontSize', 14, 'FontWeight', 'bold');

subplot(3,1,2); hold on;
plot(T, allSP_pos(:,2), 'LineWidth', 2);
plot(T, approxPos(:,2), '--', 'LineWidth', 2);
% plot(T, approxPos2(:,2), '-.', 'LineWidth', 2);
hold off; grid on;
ylabel('sp_y, non-dim');
set(gca, 'FontSize', 14, 'FontWeight', 'bold');

subplot(3,1,3); hold on;
plot(T, allSP_pos(:,3), 'LineWidth', 2);
plot(T, approxPos(:,3), '--', 'LineWidth', 2);
% plot(T, approxPos2(:,3), '-.', 'LineWidth', 2);
hold off; grid on;
ylabel('sp_z, non-dim');
set(gca, 'FontSize', 14, 'FontWeight', 'bold');

figure(); hold on;
plot((T - T0)*bcSys.charT/24/3600, (allSP_pos(:,1) - approxPos(:,1))*bcSys.charL, 'LineWidth', 2);
plot((T - T0)*bcSys.charT/24/3600, (allSP_pos(:,2) - approxPos(:,2))*bcSys.charL, 'LineWidth', 2);
plot((T - T0)*bcSys.charT/24/3600, (allSP_pos(:,3) - approxPos(:,3))*bcSys.charL, 'LineWidth', 2);
hold off; grid on;
xlabel('Days Since Base Epoch');
ylabel('Approx. Distance Error, km');
legend('x', 'y', 'z');
set(gca, 'FontSize', 14, 'FontWeight', 'bold');