% Check out corrected LPF orbit from C++
clear; clc; close all;
colors = parula(3);
bcSys = bcr4bpR_getSysParam('sun', 'earth', 'moon');

% n = load('Corrected_LPF_Sample.mat');
n = load('MinDV_LPF_Sample.mat');
o = bcr4bpR_MS2OrbStruct(n);

h_fig = open('LPF_Traj_Working_Uncorrected.fig');
figChild = allchild(h_fig);
h_ax = figChild(end);
h_uncor = h_ax.Children(end);

figure(h_fig); hold on;
h_cor = plot3(o.State(:,1), o.State(:,2), o.State(:,3), 'Color', colors(2,:), 'LineWidth', 1.5);
plot3(n.Nodes(:,1), o.Nodes(:,2), o.Nodes(:,3), 'd', 'Color', colors(2,:));
hold off; axis auto;
legend([h_uncor, h_cor], 'Uncorrected', 'Corrected');
title('LPF Traj with SP Range Constraint');

% Compute delta V on corrected orbit
totalDV = 0;
for i = 1:(length(n.Epochs)-1)
    arc = bcr4bpR_simNL(n.Nodes(i,:), 'sun', 'earth', 'moon', 't0', n.Epochs(i),...
        'tf', n.TOFs(i), 'Theta0', o.Theta0, 'Phi0', o.Phi0, 'Gamma', o.Gamma);
    DV = norm(arc.State(end,4:6) - n.Nodes(i+1,4:6));
    fprintf('DV at Node %d = %.5f m/s\n', i+1, DV*1000*bcSys.charL/bcSys.charT);
    totalDV = totalDV + DV;
end
fprintf('-----------------\nTotal DV = %.5f m/s\n', totalDV*1000*bcSys.charL/bcSys.charT);

% Check SP distance
spPos1 = bcr4bpR_locateSP(n.Epochs(10), 0, n.Theta0, n.Phi0, n.Gamma, bcSys);
spPos2 = bcr4bpR_locateSP(n.Epochs(end), 0, n.Theta0, n.Phi0, n.Gamma, bcSys);

dist1 = norm(spPos1.' - n.Nodes(10,1:3));
dist2 = norm(spPos2.' - n.Nodes(end,1:3));

fprintf('SP Distance at Node 10 is %.4e km\n', dist1*bcSys.charL);
fprintf('SP Distance at Node end is %.4e km\n', dist2*bcSys.charL);