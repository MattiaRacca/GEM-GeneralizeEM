%%  GEM - toy examples
%   Author: Mattia Racca
%   Date:   27/11/2015
clear all; close all; clc;
set(0,'DefaultFigureWindowStyle','docked');
rng(13245);

experiment_number = 6;

switch experiment_number
    case 1
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
range_points = linspace(-10 , 15)';
PDF = zeros(length(range_points),length(distribution_type));
sum_PDF = zeros(length(range_points),1);

for k=1:length(distribution_type)
    PDF(:,k) = distribution{k}.prior(end)*Gaussian(distribution{k}.mu(end), distribution{k}.sigma(end), range_points);
end
sum_PDF = computeProbabilityMixture( distribution, range_points, iterations );

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

% disp('Likelihood mixture ex.1');
% disp(computeLikelihood(distribution, data, 1));
pause;

    case 2
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
range_points = linspace(0 , 30)';
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

% disp('Likelihood mixture ex.2');
% disp(computeLikelihood(distribution, data, 1));
pause;

    case 3
%%%%%%%%%%%%%%%%% Example 3: Mixture (1 Exponential + 2 Normal) vs Mixture
%%%%%%%%%%%%%%%%% (1 Exponential + 1 Normal) vs Mixture (4 Normal)

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

%%  Mixture fitting
iterations = 100;
distribution_type = 2;

[ distribution2 ] = gem( data, iterations, distribution_type);

%%  Mixture fitting
iterations = 100;
distribution_type = [2 2 2 2]';

[ distribution3 ] = gem( data, iterations, distribution_type);

%%  Plots
range_points = linspace(0 , 30)';

% First mixture
PDF = zeros(length(range_points),length(distribution));
sum_PDF = zeros(length(range_points),1);

for k=1:length(distribution)
    switch distribution{k}.type
        case 1
            PDF(:,k) = distribution{k}.prior(end)*Exponential(distribution{k}.mu(end), range_points);
        case 2
            PDF(:,k) = distribution{k}.prior(end)*Gaussian(distribution{k}.mu(end), distribution{k}.sigma(end), range_points);
    end
    sum_PDF = sum_PDF + PDF(:,k);
end

% Second mixture
PDF2 = zeros(length(range_points),length(distribution2));
sum_PDF2 = zeros(length(range_points),1);

for k=1:length(distribution2)
    switch distribution2{k}.type
        case 1
            PDF2(:,k) = distribution2{k}.prior(end)*Exponential(distribution2{k}.mu(end), range_points);
        case 2
            PDF2(:,k) = distribution2{k}.prior(end)*Gaussian(distribution2{k}.mu(end), distribution2{k}.sigma(end), range_points);
    end
    sum_PDF2 = sum_PDF2 + PDF2(:,k);
end

% Third mixture
PDF3 = zeros(length(range_points),length(distribution3));
sum_PDF3 = zeros(length(range_points),1);

for k=1:length(distribution3)
    switch distribution3{k}.type
        case 1
            PDF3(:,k) = distribution3{k}.prior(end)*Exponential(distribution3{k}.mu(end), range_points);
        case 2
            PDF3(:,k) = distribution3{k}.prior(end)*Gaussian(distribution3{k}.mu(end), distribution3{k}.sigma(end), range_points);
    end
    sum_PDF3 = sum_PDF3 + PDF3(:,k);
end

figure
hold on;
title('Distribution and Mixture modelling: example 3')

h = histogram(data,'Normalization','pdf');
h.NumBins = 100;

pause;

for k=1:length(distribution)
    plot(range_points, PDF(:,k), 'b--', 'linewidth', 1);
end
for k=1:length(distribution2)
    plot(range_points, PDF2(:,k), 'r--', 'linewidth', 1);
end
for k=1:length(distribution3)
    plot(range_points, PDF3(:,k), 'g--', 'linewidth', 1);
end

pause;

plot(range_points, sum_PDF, 'b', 'linewidth', 2);
plot(range_points, sum_PDF2, 'r', 'linewidth', 2);
plot(range_points, sum_PDF3, 'g', 'linewidth', 2);

%%  Likelihood computation

disp('LogLikelihood analysis - 3° example');
disp('LogLikelihood 1xExp + 2xGauss:')
disp(computeLikelihood(distribution, data, 1));
disp('LogLikelihood 1xGauss:')
disp(computeLikelihood(distribution2, data, 1));
disp('LogLikelihood 4xGauss:')
disp(computeLikelihood(distribution3, data, 1));

pause;

    case 4
%%%%%%%%%%%%%%%%% Example 4: BIC test

%%  Load data
x1 = exprnd(0.5*ones(4000,1));
x2 = exprnd(1*ones(200,1));
x3 = abs(4*ones(500,1) + 8.*randn(500,1));

data = [x1; x2; x3];

%%  Mixture fitting
iterations = 100;
distribution_type = [1 1]';

[ distribution1 ] = gem( data, iterations, distribution_type);

%%  Mixture fitting
iterations = 100;
distribution_type = [1 2]';

[ distribution2 ] = gem( data, iterations, distribution_type);

%%  Plots
range_points = linspace(0 , 30)';

% First mixture
PDF1 = zeros(length(range_points),length(distribution1));

for k=1:length(distribution1)
    switch distribution1{k}.type
        case 1
            PDF1(:,k) = distribution1{k}.prior(end)*Exponential(distribution1{k}.mu(end), range_points);
        case 2
            PDF1(:,k) = distribution1{k}.prior(end)*Gaussian(distribution1{k}.mu(end), distribution1{k}.sigma(end), range_points);
    end
end
sum_PDF1 = computeProbabilityMixture( distribution1, range_points, iterations);

% Second mixture
PDF2 = zeros(length(range_points),length(distribution2));
sum_PDF2 = zeros(length(range_points),1);

for k=1:length(distribution2)
    switch distribution2{k}.type
        case 1
            PDF2(:,k) = distribution2{k}.prior(end)*Exponential(distribution2{k}.mu(end), range_points);
        case 2
            PDF2(:,k) = distribution2{k}.prior(end)*Gaussian(distribution2{k}.mu(end), distribution2{k}.sigma(end), range_points);
    end
    sum_PDF2 = sum_PDF2 + PDF2(:,k);
end

figure
hold on;
title('Distribution and Mixture modelling: example 4')

h = histogram(data,'Normalization','pdf');

pause;

for k=1:length(distribution1)
    plot(range_points, PDF1(:,k), 'b--', 'linewidth', 1);
end
for k=1:length(distribution2)
    plot(range_points, PDF2(:,k), 'r--', 'linewidth', 1);
end

pause;

plot(range_points, sum_PDF1, 'b', 'linewidth', 2);
plot(range_points, sum_PDF2, 'r', 'linewidth', 2);

%% BIC

disp('LogLikelihood analysis - 4° example');
disp('BIC 2xExp:')
disp(computeBIC(distribution1, data));
disp('BIC 1xExp + 1xGauss:')
disp(computeBIC(distribution2, data));

    case 5
%%%%%%%%%%%%%%%%% Example 5: Laplace test

%%  Load data
x1 = exprnd(0.5*ones(4000,1));
x2 = abs(- exprnd(1*ones(800,1)) + 4*ones(800,1));
x3 = exprnd(1*ones(800,1)) + 4*ones(800,1);

data = [x1; x2; x3];

%%  Mixture fitting
iterations = 100;
distribution_type = [1 3]';

[ distribution1 ] = gem( data, iterations, distribution_type);

iterations = 100;
distribution_type = [1 2]';

[ distribution2 ] = gem( data, iterations, distribution_type);

%%  Plots
range_points = linspace(0 , 10)';

% First mixture
PDF1 = zeros(length(range_points),length(distribution1));

for k=1:length(distribution1)
    PDF1(:,k) = computeProbability( distribution1, k, range_points);
end
sum_PDF1 = computeProbabilityMixture( distribution1, range_points, iterations);

% Second mixture
PDF2 = zeros(length(range_points),length(distribution2));

for k=1:length(distribution2)
    PDF2(:,k) = computeProbability( distribution2, k, range_points);
end
sum_PDF2 = computeProbabilityMixture( distribution2, range_points, iterations);


figure
hold on;
title('Distribution and Mixture modelling: example 5')

h = histogram(data,'Normalization','pdf');

pause;

for k=1:length(distribution1)
    plot(range_points, PDF1(:,k), 'b--', 'linewidth', 1);
end
for k=1:length(distribution2)
    plot(range_points, PDF2(:,k), 'r--', 'linewidth', 1);
end

pause;

plot(range_points, sum_PDF1, 'b', 'linewidth', 2);
plot(range_points, sum_PDF2, 'r', 'linewidth', 2);

    case 6
%%%%%%%%%%%%%%%%% Example 5: Laplace test

%%  Load data
x1 = trnd(123, 800,1) + 10*ones(800,1);
x2 = trnd(4, 800,1) - 20*ones(800,1);
data = [x1; x2];

%%  Mixture fitting
iterations = 100;
distribution_type = [4 4]';

initialization{1}.mu = 7;
initialization{1}.lambda = 1;
initialization{1}.nu = 10;

initialization{2}.mu = -14;
initialization{2}.lambda = 1;
initialization{2}.nu = 7;

init_vect = [1 1];

[ distribution1 ] = gem( data, iterations, distribution_type, initialization, init_vect);

%%  Plots
range_points = linspace(0 , 30)';

% First mixture
PDF1 = zeros(length(range_points),length(distribution1));

for k=1:length(distribution1)
    PDF1(:,k) = computeProbability( distribution1, k, range_points);
end
sum_PDF1 = computeProbabilityMixture( distribution1, range_points, iterations);

figure
hold on;
title('Student-t modelling: example 6')

h = histogram(data,'Normalization','pdf');

pause;

for k=1:length(distribution1)
    plot(range_points, PDF1(:,k), 'b--', 'linewidth', 1);
end

pause;

plot(range_points, sum_PDF1, 'b', 'linewidth', 2);
        
end
