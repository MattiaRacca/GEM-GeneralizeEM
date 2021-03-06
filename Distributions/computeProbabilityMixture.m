function [ P ] = computeProbabilityMixture( distribution, data, i )
%COMPUTEPROBABILITY Computes probability over data of the entire mixture
% (with i-th parameters)
n = length(data);
P = zeros(n,1);

if length(distribution) == 1
    P = (distribution{1}.prior(end) * computeProbability( distribution, 1, data));
else
    if nargin < 3
        for k=1:length(distribution)
            P = P + (distribution{k}.prior(end) * computeProbability( distribution, k, data));
        end
    else
        for k=1:length(distribution)
            P = P + (distribution{k}.prior(i) * computeProbability( distribution, k, data, i ));
        end
    end
end
end

