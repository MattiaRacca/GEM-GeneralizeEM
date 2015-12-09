function [ mu, lambda, nu] = fitStudent( data, W )
%FITSTUDENT Maximum likelihood estimate for Student-t distribution
%   Use W for setting weights to data

if nargin == 1
    [mu, temp] = fitGaussian(data);
    lambda = 1/temp;
    nu = 10;
    start = 1;
    
    for i=1:100
        u = ((nu+1)*ones(size(data))./(nu+lambda*(data-mu).^2));
        mu_new = u'*data/sum(u);
        lambda_new = mean(u.*(data-mu_new).^2)^-1;
        
        f = @(nu) -psi(nu/2) + log(nu) + 1 + psi((nu+1)/2) - log(nu+1) + mean(log(u)-u);
        if start
            start = 0;
            nu_new = fzero(f, [1E-4,1E+4]);
        else
            factor = 2;
            minmax = [nu/2, nu*2];
            while f(minmax(1))*f(minmax(2)) >= 0
                minmax(1) = minmax(1) / factor;
                minmax(2) = minmax(2) * factor;
            end
            nu_new = fzero(f, minmax);
        end
        nu = nu_new;
        mu = mu_new;
        lambda = lambda_new;
    end
else
    [mu, temp] = fitGaussian( data, W );
    lambda = 1/temp;
    nu = 0.1;
    start = 1;
    data = data.*W;
    
    for i=1:100
        u = (nu+1)/(nu+lambda*(data-mu)^2);
        mu_new = sum(data.*u)/sum(u);
        lambda_new = 1/mean(u * (data-mu_new)^2);
        
        f = @(nu) -psi(nu/2) + log(nu) + 1+psi((nu+1)/2)-log(nu+1)+mean(log(u)-u);
        if start
            start = 0;
            nu_new = fzero(f, [1E-4,1E+4]);
        else
            factor = 2;
            minmax = [nu/2, nu*2];
            while f(minmax(0))*f(minmax(1)) >= 0
                minmax(0) = minmax(0) / factor;
                minmax(1) = minmax(1) / factor;
                nu_new = fzero(f, [minmax(0),minmax(1)]);
            end
        end
        nu = nu_new;
        mu = mu_new;
        lambda = lambda_new;
    end
    
    mu = mu * length(data)/sum(W);
    lambda = lambda * sum(W)/length(data);
end
end

