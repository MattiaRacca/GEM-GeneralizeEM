function [ y ] = Exponential( mu, x )
%EXPONENTIAL Exponential distribution
%   http://se.mathworks.com/help/stats/exponential-distribution.html
y = (1/mu)*exp(-x/mu);
y = pdf('Exponential',x,mu);
end

