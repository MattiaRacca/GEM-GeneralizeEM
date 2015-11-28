function [ P ] = computeProbability( distribution, k, i, data )
%COMPUTEPROBABILITY Computes probability over data using K distribution
% (with i-th parameters) from the mixture
n = length(data);
P = zeros(n,1);

switch distribution{k}.type
    case 1
        P = Exponential(distribution{k}.mu(i), data);                          % Probability Exponential
    case 2
        P = Gaussian(distribution{k}.mu(i),distribution{k}.sigma(i), data);  % Probability Gaussian
end
end

