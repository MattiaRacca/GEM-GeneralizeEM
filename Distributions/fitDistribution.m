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
        case 3
            [fitted_distribution{k}.mu(i), fitted_distribution{k}.sigma(i)] = fitLaplace(data);
        case 4
            [fitted_distribution{k}.mu(i), fitted_distribution{k}.lambda(i), fitted_distribution{k}.nu(i)] = fitStudent(data, 0.1);
    end
else
    switch distribution{k}.type
        case 1
            fitted_distribution{k}.mu(i) = fitExponential(data,W);
        case 2
            [fitted_distribution{k}.mu(i), fitted_distribution{k}.sigma(i)] = fitGaussian(data,W);
        case 3
            [fitted_distribution{k}.mu(i), fitted_distribution{k}.sigma(i)] = fitLaplace(data,W);
        case 4
            [fitted_distribution{k}.mu(i), fitted_distribution{k}.lambda(i), fitted_distribution{k}.nu(i)] = fitStudent(data,fitted_distribution{k}.nu(i-1), W, fitted_distribution{k}.mu(i-1),fitted_distribution{k}.lambda(i-1));
    end
end

end

