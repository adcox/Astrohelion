%% Plot all the arcs from the LPF manifold propagation script

clear; clc; close all;
defineConstants;
colors = lines(5);
fontSize = 18;
lineWidth = 2;
bcSys = bcr4bpR_getSysParam('sun', 'earth', 'moon');

%   %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%   Promising Candidates from Halo Orbit
%   %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% dataDir = 'LPF_4B_NaturalManifolds/';
%
%   9, 11, 12, 57, 58, 67, 78, 84       - Nails a lunar flyby to capture in the EM system
%   51-54       - Conventional transfers, all equally promising
%   59-65, 73, 80-83, 85, 87  - Loops (large) in and out of EM system; may
%       be easier to correct
%   Other arcs
%   21-50   - Nice ribbon of manifolds that flyby and depart
%
%   %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%   Promising Candidates from Quasi-Halo Orbit #1 (Az = 150,000 km)
%   %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
dataDir = 'LPF_QH_4B_NaturalManifolds/';
%
%   Most promising candidates: 29, 34, 35, 37
%   26-35, 37-40: two pass
%   36, 83: Tightly captured
%   38: Loosely captured
%   37, 45, 51, 56, 57, 61, 84: multi-pass, perhaps easy to tweak
%   1-25: Depart through L1
%
%   %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%   Promising Candidates from Quasi-Halo Orbit #2 (Az = 500,000 km)
%   %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% dataDir = 'LPF_QH2_4B_NaturalManifolds/';
%
%   87,88,90 Very, very loosely captured
%   25-34, 43-50, 55-56, 58-59 Loose capture
%   40-42 Loose capture with large L1/L2 loops
%   51-54 Tight capture
%
%   Other arcs
%   1-8, 10-24, 89 depart through L1
%   35-39, 57, 60-62, 74-86 depart through L2
%   9 converged to strange solution
%   63-73 crash
%
%   %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%   Promising Candidates from Quasi-Halo Orbit #3 (Az = 500,000 km, 400
%   manifolds)
%   %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% dataDir = 'LPF_QH3_4B_NaturalManifolds/';
%
%   229 - like the one from Kemble's paper
%   256 - near heteroclinic connection
%
%   99-140, 151-200, 221-228: Loosely captured
%   226-246: variety of loosely captured geometries
%   141-150: Two pass and depart through L2, but too far out-of-plane
%   201-220: Tightly captured (211, 212, 217 all have one close pass)
%   
%   Other arcs
%   1-98: exit through L1
%   247-345: exit through L2
%   %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

files = dir([dataDir, 'Traj*_SEM.mat']);
ix0 = 1; %1000; % start plotting data partway through to avoid some of the halo

h_fig = figure(); hold on;
% for i = 1:length(files)
for i = [37, 56, 61]
    o = load([dataDir, files(i).name]);
    % Plot all same color
    plot3(o.State(ix0:end,1)*bcSys.charL, o.State(ix0:end,2)*bcSys.charL,...
        o.State(ix0:end,3)*bcSys.charL, 'Color', colors(1,:), 'LineWidth', lineWidth);

    % Plot all different colors
%     plot3(o.State(ix0:end,1)*bcSys.charL, o.State(ix0:end,2)*bcSys.charL,...
%         o.State(ix0:end,3)*bcSys.charL, 'LineWidth', lineWidth);
    
    axis equal;
end
grid on; limits = axis;
bcr4bpR_plotEnv('Sun', 'Earth', 'Moon', 't0', o.Time(1), 'tf', o.Time(end),...
    'scale', bcSys.charL, 'plotsp', 'on', 'LineWidth', lineWidth, 'surf', 'off');
hold off; axis(limits);
xlabel('x, km');
ylabel('y, km');
zlabel('z, km');
set(gca, 'FontSize', fontSize, 'FontWeight', fontWeight);



%% Detailed Notes for Halo Manifolds
%
%   9 - Captures VERY well but has significant out-of-plane components. If
%       I can control the flyby, I might be able to tame the out-of-plane bits. 
%   67 - Cool flow pattern that loops around the 3D lunar orbit; probably
%       too out of plane for good SP encounters
%
%   81, 82, 83, 87 - Has very close passes - COULD BE VERY GOOD CANDIDATE
%   84 - Several flybys look close to the right in-plane amplitude

%% Detailed Notes for Quasi-Halo Manifolds
%
%   34 is the most similar to the feasible solution Kemble found. Try
%   correcting this one to have one SP encounter at the last opportunitiy;
%   the design Kemble shows loops back for one more encounter
%
%   Manifolds between 36 and 37 ought to loop around L2 and return; 36
%   loops before reaching L2, and 37 *almost* turns back, but doesn't quite
