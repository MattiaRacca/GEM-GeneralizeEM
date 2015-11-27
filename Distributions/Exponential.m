function [ y ] = Exponential( mu, x )
%EXPONENTIAL Exponential distribution
y = (1/mu)*exp(-x/mu);
end

