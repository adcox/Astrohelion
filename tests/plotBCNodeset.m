clear; clc;
o = load('CorrectedSet.mat');
% o = load('HalfLyap.mat');

figure(); hold on;
totalDV = 0;
arcs = [];
for i = 1:(length(o.TOFs)-1)
    arc = bcr4bpR_simNL(o.Nodes(i,:), o.P1, o.P2, o.P3, ...
        't0', o.Epochs(i), 'tf', o.TOFs(i));
    totalDV = totalDV + norm(o.Nodes(i+1,4:6) - arc.State(end,4:6));
    plot3(arc.State(:,1), arc.State(:,2), arc.State(:,3), 'LineWidth', 2);
    arcs = [arcs; arc];
end
plot3(o.Nodes(end,1), o.Nodes(end,2), o.Nodes(end,3), 'sk');
axis equal; limits = axis;
bcr4bpR_plotEnv(o.P1, o.P2, o.P3, o.Epochs(1), o.Epochs(end));
hold off; grid on; axis(limits);
fprintf('Total Delta-V = %.4e\n', totalDV);