function [ distribution ] = gem( data, iterations, distribution_type, initialization)
%GEM General Expectation Maximization algorithm
%   [DISTRIBUTION] = GEM(DATA, DISTRIBUTION_TYPE, ITERATION) fits mixture of different
%   distributions to 1-D data.  
%   DISTRIBUTION_TYPE specifies the mixture nature, i.e. the number of
%   distributions and their type. The i-th number in the vector is the type
%   of the distribution.
%       e.g.) [1 2 2] fits one exponential and two normal distributions
%   
%   Supported distributions:
%       1) Exponential          http://se.mathworks.com/help/stats/exponential-distribution.html
%       2) Normal (Gaussian)    http://se.mathworks.com/help/stats/normal-distribution.html
%
%   Author: Mattia Racca        Date: 27/11/2015

addpath('Distributions');

N = length(data);               % Number of samples
K = length(distribution_type);  % Number of components in the mixture

% Integrate information about type in distribution
for k=1:K
    distribution{k}.type = distribution_type(k);
end

if K == 1
    switch distribution_type
        case 1
            distribution{K}.mu = fitExponential(data);
        case 2
            [distribution{K}.mu, distribution{K}.sigma] = fitGaussian(data);
    end
else
    %% Initialization of the distributions and priors
    
    % Prior initialization
    for k=1:K
        distribution{k}.prior(1) = 1/K;
    end
    
    % Parameters initialization
    if nargin >= 4
        %   initialization from proper parameter
        for k=1:K
            switch distribution_type(k)
                case 1
                    distribution{k}.mu(1) = initialization{k}.mu;
                case 2
                    distribution{k}.mu(1) = initialization{k}.mu;
                    distribution{k}.sigma(1) = initialization{k}.sigma;
            end
        end
    else
        %   random initialization
        for k=1:K
            switch distribution_type(k)
                case 1
                   distribution{k}.mu(1) = fitExponential(randsample(data,floor(length(data)/K)));
                case 2
                   [distribution{k}.mu(1), distribution{k}.sigma(1)] = fitGaussian(randsample(data,floor(length(data)/K)));
            end
        end
    end
    
    %% EM algorithm
    
    P    = zeros(N,K);                 % Probabilities
    Post = zeros(N,K);                 % Posteriors
    
    for i=2:iterations
        % Compute probabilities
        for k=1:K
            switch distribution_type(k)
                case 1
                    P(:,k) = Exponential(distribution{k}.mu(i-1), data);                          % Probability Exponential
                case 2
                    P(:,k) = Gaussian(distribution{k}.mu(i-1),distribution{k}.sigma(i-1), data);  % Probability Gaussian
            end
        end
        
        % Compute posteriors
        sum_posteriors = zeros(N,1);
        for k=1:K
            Post(:,k) = P(:,k)*distribution{k}.prior(i-1);
            sum_posteriors = sum_posteriors + Post(:,k);
        end
        
        for k=1:K
            Post(:,k) = Post(:,k) ./ sum_posteriors;
        end
        
        % Reestimation of parameters
        for k=1:K
            switch distribution_type(k)
                case 1
                    distribution{k}.mu(i) = fitExponential(data, Post(:,k));
                case 2
                    [distribution{k}.mu(i), distribution{k}.sigma(i)] = fitGaussian(data, Post(:,k));
            end
        end
        
        % Reestimation of priors
        for k=1:K
            distribution{k}.prior(i) = mean(Post(:,k));
        end
    end
end
end