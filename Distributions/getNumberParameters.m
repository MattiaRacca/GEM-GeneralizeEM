function [ K ] = getNumberParameters( distribution )
%GETNUMBERPARAMETERS Gives the number of parameters for BIC computation
K = 0;
switch distribution.type
        case 1
            K = K + 1;
        case 2
            K = K + 2;
end
end

