function [ mean_exp ] = fitExponential( data, W )
%FITEXPONENTIAL Maximum likelihood estimate for Exponential
%   Use W for setting weights to data

if nargin == 1
    mean_exp = mean(data);
else
    mean_exp = W'*data / sum(W);
end

end

