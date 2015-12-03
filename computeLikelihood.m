function [ likelihood ] = computeLikelihood( distribution, data, log_likelihood)
%COMPUTELIKELIHOOD Computes (log)likelihood of a mixture given data
%   [L] = COMPUTELIKELIHOOD(M, D, log) computes L likelihood, given the
%   mixture M and the data D. If log == true, it computes the loglikelihood

N = length(data);               % Number of samples
K = length(distribution);       % Number of components in the mixture

if(log_likelihood)
    likelihood = 0;                 % LogLikelihood
    for i=1:N
        sample_L = 0;
        for k=1:K
            sample_L = sample_L + distribution{k}.prior(end)*computeProbability( distribution, k, data(i));
        end
        likelihood = likelihood + log(sample_L);
    end
else
    likelihood = 1;                 % Likelihood
    for i=1:N
        sample_L = 0;
        for k=1:K
            sample_L = sample_L + distribution{k}.prior(end)*computeProbability( distribution, k, data(i));
        end
        likelihood = likelihood * sample_L;
    end
end


end

