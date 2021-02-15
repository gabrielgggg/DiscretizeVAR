clear;

% Solution
midIx = h5read('myMCvar.h5', '/DVAR/midIx') + 1;
grids = h5read('myMCvar.h5', '/DVAR/grids');
transitions = h5read('myMCvar.h5', '/DVAR/transitions');
varIntercept = h5read('myMCvar.h5', '/DVAR/varIntercept');
varRho = h5read('myMCvar.h5', '/DVAR/varRho');
varSigma = h5read('myMCvar.h5', '/DVAR/varSigma');

sz = size(grids, 2);

ttmp = transitions^5000;
statdist = ttmp(1, :);

T = 5000;
path = simulateMarkov(T, transitions);
path1 = grids(path, 1);
path2 = grids(path, 2);
path3 = grids(path, 3);

preciseEps = mvnrnd(zeros([sz, 1]), varSigma, T);
precisePath = zeros([T, sz]);
precisePath(1, :) = grids(1, :);
for tIx = 2:T
  precisePath(tIx, :) = (varIntercept + varRho * precisePath(tIx-1, :)' + preciseEps(tIx, :)')';
end
precise1 = precisePath(:, 1);
precise2 = precisePath(:, 2);
precise3 = precisePath(:, 3);

% scatter(grids(:, 1), grids(:, 2), statdist*1000+1);

% scatter3(grids(:, 1), grids(:, 2), grids(:, 3), rrow*1000+100);

close all;
figure;
% subplot(1, 2, 1);
% scatter3(grids(:, 1), grids(:, 2), statdist); hold on;
% scatter(path1, path2, 'rx');
% scatter3(precise1, precise2, precise3, 'rx'); hold on;
scatter3(grids(:, 1), grids(:, 2), grids(:, 3), 'bo', 'filled'); 
hold on;
scatter3(grids(midIx, 1), grids(midIx, 2), grids(midIx, 3), 'kd', 'filled'); 


k = boundary(grids);
trisurf(k,grids(:,1),grids(:,2),grids(:,3),'Facecolor','blue','FaceAlpha',0.15)

% subplot(1, 2, 2);
% scatter3(grids(:, 1), grids(:, 2), squeeze(transitions(5, :)));



