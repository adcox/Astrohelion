clear; clc; close all;

C = 3.144;
dirName = 'data';
P1 = 'earth';
P2 = 'moon';
mu = cr3bp_getSysParam(P1, P2);

load(sprintf('%s/JC_%03d_ICs.mat', dirName, round(C*100)));
ICs = ICs.';

files = dir(sprintf('%s/JC_%03d_Orb*.mat', dirName, round(C*100)));

%% Create a concatenated file
allReturns = [];
orbitNums = [];
for n = 1:size(files,1)
    filename = sprintf('%s/%s', dirName, files(n).name);
    load(filename);
    allReturns = [allReturns; data.'];
    orbitNums = [orbitNums; ones(size(data,1),1)*n];
end

%% Compute Constant Jacobi Curve for map
x = min(allReturns(:,1)):0.0001:max(allReturns(:,1));
c = 1;

for i = 1:length(x)
    d = abs(x(i) + mu);
    r = abs(x(i) - 1 + mu);
    % Find maximum value of x_dot for this x and C; we assume y-dot is at
    % its minimum mangitude (zero)
    x_dot = sqrt( -C + 2*(1-mu)/d + 2*mu/r + x(i)^2 );

    % only save points that are real
    if(x_dot == x_dot')
        constCPts(c,:) = [x(i), x_dot];
        c = c+1;
    else
        % Append a NaN point so that ZVCs aren't connect across empty space
        if(c > 1)
            if(constCPts(c-1,:) == constCPts(c-1,:))
                constCPts(c,:) = [NaN,NaN];
                c=c+1;
            end
        end
    end
end

%% Save Data
save(sprintf('%s/JC_%03d_MapData.mat', dirName, round(C*100)), ...
    'allReturns', 'ICs', 'C', 'P1', 'P2', 'constCPts',...
    'orbitNums');