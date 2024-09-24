% Read the data from the CSV file
data = readtable('response_time_heatplot.txt');

% Extract the columns into separate arrays
lambda = data.lambda;
mu = data.mu;
values = data.Avg_cluster;

% Create a grid for interpolation
[lambdaGrid, muGrid] = meshgrid(linspace(min(lambda), max(lambda), 100), linspace(min(mu), max(mu), 100));

% Interpolate the values
valuesGrid = griddata(lambda, mu, values, lambdaGrid, muGrid, 'cubic');

% Plot the heatmap
figure;
imagesc(lambdaGrid(1,:), muGrid(:,1), valuesGrid);
set(gca, 'YDir', 'normal'); % Correct the Y-axis direction
colorbar;
colormap(icefire); % Apply the 'jet' colormap
xlabel('\lambda');
ylabel('\mu');
title('Average cluster');
caxis([0 20]);
