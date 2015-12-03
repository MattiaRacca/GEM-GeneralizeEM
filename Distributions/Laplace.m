function [ y ] = Laplace(mu, sigma, x)
%LAPLACE laplace probability distribution
%   https://en.wikipedia.org/wiki/Laplace_distribution
y = exp(-abs(x-mu)/(sigma))/(sigma*2);
end

