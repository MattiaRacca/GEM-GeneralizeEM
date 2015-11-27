function [fitted_distribution] = fitDistribution( distribution, k, i, data, W)
%FITDISTRIBUTION Generic version of the various fit
%   FITDISTRO = FITDISTRIBUTION(DISTRO, K, I, DATA, W)
%   DISTRO = Mixture to be fitted
%   K = component inside the mixture to fit
%   I = FITTING STEP (inside EM algorithm)
%   DATA = data to fit
%   W = weight for data

fitted_distribution = distribution;
if nargin < 5
    switch distribution{k}.type
        case 1
            fitted_distribution{k}.mu(i) = fitExponential(data);
        case 2
            [fitted_distribution{k}.mu(i), fitted_distribution{k}.sigma(i)] = fitGaussian(data);
    end
else
    switch distribution{k}.type
        case 1
            fitted_distribution{k}.mu(i) = fitExponential(data,W);
        case 2
            [fitted_distribution{k}.mu(i), fitted_distribution{k}.sigma(i)] = fitGaussian(data,W);
    end
end

end

