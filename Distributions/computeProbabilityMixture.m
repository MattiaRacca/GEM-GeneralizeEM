function [ P ] = computeProbabilityMixture( distribution, i, data )
%COMPUTEPROBABILITY Computes probability over data of the entire mixture
% (with i-th parameters)
n = length(data);
P = zeros(n,1);

for k=1:length(distribution)
    P = P + (distribution{k}.prior(i) * computeProbability( distribution, k, i, data ))';
end

end

