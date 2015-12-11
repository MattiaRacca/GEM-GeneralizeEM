function [ K ] = getNumberParameters( distribution )
%GETNUMBERPARAMETERS Gives the number of parameters for BIC computation
K = 0;
switch distribution.type
    case 1              % Exponential
        K = K + 1;
    case {2, 3}         % Guassian and Laplace
        K = K + 2;
    case 4
        K = K + 3;      % Student-t
end
end

