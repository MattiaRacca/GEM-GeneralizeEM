function [ y ] = Gaussian(mu, sigma, x)
%GUASSIAN guassian probability distribution
%	http://se.mathworks.com/help/stats/normal-distribution.html
% y = exp(-(x-mu).^2/(2*sigma.^2))/(sigma*sqrt(2*pi));
y = pdf('Normal',x,mu,sigma);
end

