% Plot the stability index and periapse altitude
clear; close all;
defineConstants;

numOrbits = 50;

load('data/EML2_NButterfly_PeriStability');
famData = load('../share/families_natParam_checked/EM_L2_NButterfly.mat');
% load('data/EML1_Halo_PeriStability');
% famData = load('../share/families_pac_checked/EM_L1_Halo.mat');

% Make all indices positive
StabilityIx = abs(StabilityIx);

h_periVsStable = figure();
colormap(cool);
scatter(PeriAlt, StabilityIx, 25, StabilityIx, 'filled');
grid on;
xlabel('Lunar Periapsis Altitude, km', 'FontSize', fontSize, 'FontWeight', fontWeight);
ylabel('Stability Index', 'FontSize', fontSize, 'FontWeight', fontWeight);
title('Butterfly Family: Stability vs Perilune Altitude',...
    'FontSize', fontSize, 'FontWeight', fontWeight);

colors = cool(200);
colors = [colors;0 0 0];    %% Append black
orbitStep = floor(max(size(StabilityIx))/numOrbits);
mirror = eye(6);
mirror(3,3) = -1;
mirror(6,6) = -1;

h_family_stability = figure(); hold on;
colormap(cool);
for i = 1:orbitStep:max(size(StabilityIx))
    if(~isnan(StabilityIx(i)))
        col_ix = floor((StabilityIx(i) - min(StabilityIx))/range(StabilityIx)*(max(size(colors)) - 1));
    else
        col_ix = max(size(colors));
    end
    if(col_ix == 0), col_ix = 1; end
    fprintf('Orbit %03d: Color Index = %03d\n', i, col_ix);
    o = cr3bp_simNL(famData.MemberData(i,1:6)*mirror, famData.P1, famData.P2, ...
        'tf', famData.MemberData(i,7), 'exitCond', 'none');
    plot3(o.State(:,1), o.State(:,2), o.State(:,3), 'color', colors(col_ix,:),...
        'LineWidth', lineWidth);
end
limits = axis;
caxis([min(StabilityIx), max(StabilityIx)]);
h_c = colorbar;
cr3bp_plotEnv(famData.P1, famData.P2);
axis(limits); axis equal; grid on;
xlabel('x, non-dim', 'FontSize', fontSize, 'FontWeight', fontWeight);
ylabel('y, non-dim', 'FontSize', fontSize, 'FontWeight', fontWeight);
ylabel(h_c, 'Stability Index', 'FontSize', fontSize, 'FontWeight', fontWeight);
zlabel('z, non-dim', 'FontSize', fontSize, 'FontWeight', fontWeight);