clear;

% Solution
grids = h5read('stochVol.h5', '/DsvVAR/grids');
transitions = h5read('stochVol.h5', '/DsvVAR/transitions');
varIntercept = h5read('stochVol.h5', '/DsvVAR/varIntercept');
varRho = h5read('stochVol.h5', '/DsvVAR/varRho');
varSigma = h5read('stochVol.h5', '/DsvVAR/varSigma');

sz = size(grids, 2);

close all;
figure;
scatter3(grids(:, 1), grids(:, 2), grids(:, 3), 'bo', 'filled'); 

