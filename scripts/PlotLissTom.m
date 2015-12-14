% Plot Liss from Tom's Method

clear; clc; close all;
defineConstants;

o = load('data/Liss_TomMethod_Traj.mat');
n = load('data/Liss_TomMethod_Nodes.mat');

[mu3B,charT,charL,~] = cr3bp_getSysParam(o.P1, o.P2);
LPts = cr3bp_getEquilibPts(mu3B,'matrix', 'on');
LPts = LPts*charL;

% Convert to dimensional units
o.State = o.State*charL;
o.State(:,4:6) = o.State(:,4:6)/charT;
n.Nodes = n.Nodes*charL;
n.Nodes(:,4:6) = n.Nodes(:,4:6)/charT;

figure(); hold on;
plot3(o.State(:,1) - LPts(1,1), o.State(:,2), o.State(:,3), 'LineWidth', lineWidth);
% plot3(n.Nodes(:,1) - LPts(1,1), n.Nodes(:,2), n.Nodes(:,3), 'ok');
% for i = 1:size(n.Nodes,1)
%     text(n.Nodes(i,1) - LPts(1,1), n.Nodes(i,2), n.Nodes(i,3), int2str(i));
% end
hold off; axis equal; grid on;
xlabel('x');
ylabel('y');
zlabel('z');
set(gca, 'FontSize', fontSize, 'FontWeight', fontWeight);