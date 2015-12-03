function [ mean_lap, std_lap ] = fitLaplace( data, W )
%FITLAPLACE fits laplace distribution based on data
%   Use W for setting weights to data

if nargin == 1
    mean_lap = median(data);
    std_lap  = mean(abs(data - mean_lap));
else
    mean_lap = W'*data / sum(W);
    std_lap  = W'*(abs(data - mean_lap))/sum(W);
end

end

