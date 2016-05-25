clear; clc;
o = load('emDRO_newNodes.mat');

figure(); hold on;
totalDV = 0;
arcs = [];
for i = 1:(length(o.TOFs)-1)
    arc = cr3bp_simNL(o.Nodes(i,:), o.P1, o.P2, 'tf', o.TOFs(i));
    totalDV = totalDV + norm(o.Nodes(i+1,4:6) - arc.State(end,4:6));
    plot(arc.State(:,1), arc.State(:,2), 'LineWidth', 2);
    arcs = [arcs; arc];
end
axis equal;
plot(o.Nodes(:,1), o.Nodes(:,2), 'kd', 'MarkerFaceColor', 'k');
limits = axis;
cr3bp_plotEnv(o.P1, o.P2);
hold off; grid on; axis(limits);
xlabel('x, nondim');
ylabel('y, nondim');
set(gca, 'FontSize', 18, 'FontWeight', 'bold');

fprintf('Total Delta-V = %.4e\n', totalDV);
