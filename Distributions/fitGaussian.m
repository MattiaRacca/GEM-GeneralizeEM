function [ mean_gauss, std_gauss ] = fitGaussian( data, W )
%FITGAUSSIAN Maximum likelihood estimate for Normal distribution
%   Use W for setting weights to data

if nargin == 1
    mean_gauss = mean(data);
    std_gauss  = std(data);
else
    mean_gauss = W'*data / sum(W);
    std_gauss  = sqrt((W'*(data-mean_gauss).^2)/sum(W));
end

end

