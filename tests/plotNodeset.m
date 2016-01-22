clear; clc;
o = load('correctedSet.mat');

figure(); hold on;
totalDV = 0;
arcs = [];
for i = 1:(length(o.TOFs)-1)
    arc = cr3bp_simNL(o.Nodes(i,:), 'earth', 'moon', 'tf', o.TOFs(i));
    totalDV = totalDV + norm(o.Nodes(i+1,4:6) - arc.State(end,4:6));
    plot(arc.State(:,1), arc.State(:,2), 'LineWidth', 2);
    arcs = [arcs; arc];
end
plot(o.Nodes(end,1), o.Nodes(end,2), 'sk');
hold off; grid on; axis equal;
fprintf('Total Delta-V = %.4e\n', totalDV);