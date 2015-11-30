function [ P ] = computeProbability( distribution, k, data, i)
%COMPUTEPROBABILITY Computes probability over data using K distribution
% (with i-th parameters) from the mixture
n = length(data);
P = zeros(n,1);

if nargin > 3
    switch distribution{k}.type
        case 1
            P = Exponential(distribution{k}.mu(i), data);                          % Probability Exponential
        case 2
            P = Gaussian(distribution{k}.mu(i),distribution{k}.sigma(i), data);  % Probability Gaussian
    end
else
    switch distribution{k}.type
        case 1
            P = Exponential(distribution{k}.mu(end), data);                          % Probability Exponential
        case 2
            P = Gaussian(distribution{k}.mu(end),distribution{k}.sigma(end), data);  % Probability Gaussian
    end
end
end

