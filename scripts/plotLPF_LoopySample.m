% Check out corrected LPF orbit from C++
clear; clc; close all;
defineConstants;
colors = parula(10);
lineWidth = 2;
fontSize = 18;

bcSys = bcr4bpR_getSysParam('sun', 'earth', 'moon');
L = bcSys.charL;

N_uncor = load('LPF_LoopySample_Uncorrected');
N_cor = load('LPF_LoopySample_minDV');

O_uncor = bcr4bpR_MS2OrbStruct(N_uncor);
O_cor = bcr4bpR_MS2OrbStruct(N_cor);

%%
figure(1); hold on;
plot(O_uncor.State(:,1)*L, O_uncor.State(:,2)*L, 'Color', colors(3,:), 'LineWidth', lineWidth);
plot(O_cor.State(:,1)*L, O_cor.State(:,2)*L, 'Color', colors(5,:), 'LineWidth', lineWidth);
axis equal;
plot(N_uncor.Nodes(:,1)*L, N_uncor.Nodes(:,2)*L, 'kd', 'MarkerSize', 7);
plot(N_cor.Nodes([6,13],1)*L, N_cor.Nodes([6,13],2)*L, 'mv', 'MarkerFaceColor', 'm', 'MarkerSize', 8);
limits = axis;
envData = bcr4bpR_plotEnv('sun', 'earth', 'moon', 't0', O_uncor.Time(1),...
    'tf', O_uncor.Time(end), 'theta0', O_uncor.Theta0, 'phi0', O_uncor.Phi0,...
    'gamma', O_uncor.Gamma, 'plotsp', 'on', 'scale', L, 'LineWidth', lineWidth);
hold off; grid on; axis(limits);
xlabel('x, km');
ylabel('y, km');
legend('Uncorrected', 'Corrected', 'Location', 'NorthWest');
set(gca, 'FontSize', fontSize, 'FontWeight', fontWeight);
h = findobj(gca, 'Type', 'Surface');
delete(h);

%%
figure(2); hold on;
plot(envData.P3Pos(:,1), envData.P3Pos(:,3), 'Color', col_gray, 'LineWidth', lineWidth);
plot(0, 0, 'ok', 'MarkerFaceColor', 'k');
plot(envData.SPPos(end,1), envData.SPPos(end,3), 'r*', 'MarkerSize', 15);
h_uncor = plot(O_uncor.State(:,1)*L, O_uncor.State(:,3)*L, 'Color', colors(3,:), 'LineWidth', lineWidth);
h_cor = plot(O_cor.State(:,1)*L, O_cor.State(:,3)*L, 'Color', colors(5,:), 'LineWidth', lineWidth);
plot(N_uncor.Nodes(:,1)*L, N_uncor.Nodes(:,3)*L, 'kd', 'MarkerSize', 7);
plot(N_cor.Nodes([6,13],1)*L, N_cor.Nodes([6,13],3)*L, 'mv', 'MarkerFaceColor', 'm', 'MarkerSize', 8);
axis equal; hold off; grid on;
xlabel('x, km');
ylabel('z, km');
legend([h_uncor, h_cor], 'Uncorrected', 'Corrected', 'Location', 'NorthWest');
set(gca, 'FontSize', fontSize, 'FontWeight', fontWeight);

%%

% Compute delta V on corrected orbit
totalDV = 0;
for i = 1:(length(N_cor.Epochs)-1)
    arc = bcr4bpR_simNL(N_cor.Nodes(i,:), 'sun', 'earth', 'moon', 't0', N_cor.Epochs(i),...
        'tf', N_cor.TOFs(i), 'Theta0', O_cor.Theta0, 'Phi0', O_cor.Phi0, 'Gamma', O_cor.Gamma);
    DV = norm(arc.State(end,4:6) - N_cor.Nodes(i+1,4:6));
    fprintf('DV at Node %d = %.5f m/s\n', i+1, DV*1000*bcSys.charL/bcSys.charT);
    totalDV = totalDV + DV;
end
fprintf('-----------------\nTotal DV = %.5f m/s\n', totalDV*1000*bcSys.charL/bcSys.charT);

% Check SP distance
spPos1 = bcr4bpR_locateSP(N_cor.Epochs(10), 0, N_cor.Theta0, N_cor.Phi0, N_cor.Gamma, bcSys);
spPos2 = bcr4bpR_locateSP(N_cor.Epochs(end), 0, N_cor.Theta0, N_cor.Phi0, N_cor.Gamma, bcSys);

dist1 = norm(spPos1.' - N_cor.Nodes(10,1:3));
dist2 = norm(spPos2.' - N_cor.Nodes(end,1:3));

fprintf('SP Distance at Node 10 is %.4e km = %.4fmm\n', dist1*bcSys.charL, dist1*bcSys.charL*1e6);
fprintf('SP Distance at Node end is %.4e km = %.4fmm\n', dist2*bcSys.charL, dist2*bcSys.charL*1e6);