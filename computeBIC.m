function [ BIC ] = computeBIC( distribution,  data, loglikelihood)
%COMPUTEBIC Compute Bayesian Information Criterion
%   https://en.wikipedia.org/wiki/Bayesian_information_criterion
%
%   BIC = COMPUTEBIC(DISTR, DATA) computes first the logL and then BIC
%   BIC = COMPUTEBIC(DISTR, DATA, LOGL) computes BIC using the passed logL

N = length(data);

%   Compute number of parameters K
K = 0;
for i=1:length(distribution)
    switch distribution{i}.type
        case 1
            K = K + 1;
        case 2
            K = K + 2;
    end
end

K = K + length(distribution) -1;

if nargin > 2
    BIC = -2*loglikelihood + K*log(N);
else
    BIC = -2*computeLikelihood(distribution, data, 1) + K*log(N);
end

end

