function [ y ] = Gaussian(mu, sigma, x)
%GUASSIAN guassian probability distribution
y = exp(-(x-mu).^2/(2*sigma.^2))/(sigma*sqrt(2*pi));
end

