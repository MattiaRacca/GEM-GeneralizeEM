%%  GEM - toy examples
%   Author: Mattia Racca
%   Date:   27/11/2015
clear all; close all; clc;
set(0,'DefaultFigureWindowStyle','docked');

%%%%%%%%%%%%%%%%% Example 1: GMM (3 Normal distribution)

%%  Load data
x1 = randn(4000,1);
x2 = 8*ones(2000,1) + randn(2000,1);
x3 = -5*ones(1000,1) + randn(1000,1);

data = [x1; x2; x3];

%%  Mixture fitting
iterations = 100;
distribution_type = [2 2 2]';

[ distribution ] = gem( data, iterations, distribution_type);

%%  Plots
range_points = linspace(-10 , 15);
PDF = zeros(length(range_points),length(distribution_type));
sum_PDF = zeros(length(range_points),1);

for k=1:length(distribution_type)
    PDF(:,k) = distribution{k}.prior(end)*Gaussian(distribution{k}.mu(end), distribution{k}.sigma(end), range_points);
    sum_PDF = sum_PDF + PDF(:,k);
end

figure
hold on;
title('Distribution and Mixture modelling: example 1')

h = histogram(data,'Normalization','pdf');
h.NumBins = 100;

pause;

for k=1:length(distribution_type)
    plot(range_points, PDF(:,k), 'r', 'linewidth', 1);
end

pause;

plot(range_points, sum_PDF, 'b', 'linewidth', 2);

pause;

%%%%%%%%%%%%%%%%% Example 2: Mixture (2 Exponential + 1 Normal)

%%  Load data
x1 = exprnd(ones(1000,1));
x2 = abs(15*ones(500,1) + randn(500,1));
x3 = exprnd(5*ones(1000,1));

data = [x1; x2; x3];

%%  Mixture fitting
iterations = 100;
distribution_type = [1 1 2]';

[ distribution ] = gem( data, iterations, distribution_type);

%%  Plots
range_points = linspace(0 , 30);
PDF = zeros(length(range_points),length(distribution_type));
sum_PDF = zeros(length(range_points),1);

for k=1:length(distribution_type)
    switch distribution{k}.type
        case 1
            PDF(:,k) = distribution{k}.prior(end)*Exponential(distribution{k}.mu(end), range_points);
        case 2
            PDF(:,k) = distribution{k}.prior(end)*Gaussian(distribution{k}.mu(end), distribution{k}.sigma(end), range_points);
    end
    sum_PDF = sum_PDF + PDF(:,k);
end

figure
hold on;
title('Distribution and Mixture modelling: example 2')

h = histogram(data,'Normalization','pdf');
h.NumBins = 100;

pause;

for k=1:length(distribution_type)
    plot(range_points, PDF(:,k), 'r', 'linewidth', 1);
end

pause;

plot(range_points, sum_PDF, 'b', 'linewidth', 2);

pause;

%%%%%%%%%%%%%%%%% Example 3: Mixture (1 Exponential + 2 Normal)

%%  Load data
x1 = exprnd(4*ones(1000,1));
x2 = abs(8*ones(500,1) + 2.*randn(500,1));
x3 = abs(4*ones(1000,1) + randn(1000,1));
x4 = abs(6*ones(300,1) + 5.*randn(300,1));

data = [x1; x2; x3; x4];

%%  Mixture fitting
iterations = 100;
distribution_type = [1 2 2]';

[ distribution ] = gem( data, iterations, distribution_type);

%%  Plots
range_points = linspace(0 , 30);
PDF = zeros(length(range_points),length(distribution_type));
sum_PDF = zeros(length(range_points),1);

for k=1:length(distribution_type)
    switch distribution{k}.type
        case 1
            PDF(:,k) = distribution{k}.prior(end)*Exponential(distribution{k}.mu(end), range_points);
        case 2
            PDF(:,k) = distribution{k}.prior(end)*Gaussian(distribution{k}.mu(end), distribution{k}.sigma(end), range_points);
    end
    sum_PDF = sum_PDF + PDF(:,k);
end

figure
hold on;
title('Distribution and Mixture modelling: example 3')

h = histogram(data,'Normalization','pdf');
h.NumBins = 100;

pause;

for k=1:length(distribution_type)
    plot(range_points, PDF(:,k), 'r', 'linewidth', 1);
end

pause;

plot(range_points, sum_PDF, 'b', 'linewidth', 2);

