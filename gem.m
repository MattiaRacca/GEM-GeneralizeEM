function [ distribution ] = gem( data, iterations, distribution_type, initialization, initializ_vect)
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
%       3) Laplace              https://en.wikipedia.org/wiki/Laplace_distribution
%       4) Student-t            http://se.mathworks.com/help/stats/tpdf.html
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
    distribution{K}.prior = 1;
    distribution = fitDistribution( distribution, K, 1, data);
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
            if(initializ_vect)
                switch distribution_type(k)
                    case 1
                        distribution{k}.mu(1) = initialization{k}.mu;
                    case 2
                        distribution{k}.mu(1) = initialization{k}.mu;
                        distribution{k}.sigma(1) = initialization{k}.sigma;
                    case 3
                        distribution{k}.mu(1) = initialization{k}.mu;
                        distribution{k}.sigma(1) = initialization{k}.sigma;
                    case 4
                        distribution{k}.mu(1) = initialization{k}.mu;
                        distribution{k}.lambda(1) = initialization{k}.lambda;
                        distribution{k}.nu(1) = initialization{k}.nu;
                end
            else
                distribution = fitDistribution( distribution, k, 1, randsample(data,floor(length(data)/K)));
            end
        end
    else
        %   random initialization
        for k=1:K
            if  distribution_type(k) == 4
                distribution{k}.mu(1) = randsample(data,1);
                distribution{k}.lambda(1) = 10*rand();
                distribution{k}.nu(1) = 0.1;
            else
                distribution = fitDistribution( distribution, k, 1, randsample(data,floor(length(data)/K)));
            end
        end
    end
    
    %% EM algorithm
    
    P    = zeros(N,K);                 % Probabilities
    Post = zeros(N,K);                 % Posteriors
    
    for i=2:iterations
        % Compute probabilities
        for k=1:K
            P(:,k) = computeProbability( distribution, k, data, i-1 );
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
            distribution = fitDistribution( distribution, k, i, data, Post(:,k));
        end
        
        % Reestimation of priors
        for k=1:K
            distribution{k}.prior(i) = mean(Post(:,k));
        end
    end
end
end