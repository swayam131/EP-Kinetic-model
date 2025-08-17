% Reading data
data = readtable('response_time_heatplot.txt');

% Extract data
lambda = data.lambda;
mu = data.mu;
values = data.Avg_cluster;

% interpolation
[lambdaGrid, muGrid] = meshgrid(linspace(min(lambda), max(lambda), 100), linspace(min(mu), max(mu), 100));


valuesGrid = griddata(lambda, mu, values, lambdaGrid, muGrid, 'cubic');

% Plot
figure;
imagesc(lambdaGrid(1,:), muGrid(:,1), valuesGrid);
set(gca, 'YDir', 'normal'); % Correct the Y-axis direction
colorbar;
colormap(icefire); % Apply the 'jet' colormap
xlabel('\lambda');
ylabel('\mu');
title('Average cluster');
caxis([0 20]);

