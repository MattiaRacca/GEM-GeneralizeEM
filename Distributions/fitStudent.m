function [ mu, lambda, nu] = fitStudent( data, nu_init, W, mu_init, lambda_init)
%FITSTUDENT Maximum likelihood estimate for Student-t distribution
%   Use W for setting weights to data

if nargin == 2
    [mu, temp] = fitGaussian(data);
    lambda = 1/temp;
    nu = nu_init;
    start = 1;
    
    for i=1:100
        
        %   A derivation of the Em Updates for Finding the ML param
        %   estimation of the Student's t Distribution
        %   C. Scheffler, 2008
        %   http://www.inference.phy.cam.ac.uk/cs482/publications/scheffler2008derivation.pdf
        
        
        u = ((nu+1)*ones(size(data))./(nu+lambda*(data-mu).^2));
        mu_new = u'*data/sum(u);
        lambda_new = mean(u.*(data-mu_new).^2)^-1;
        
        f = @(nu) -psi(nu/2) + log(nu) + 1 + psi((nu+1)/2) - log(nu+1) + mean(log(u)-u);
        if start
            start = 0;
            if(prod(isfinite([f(1E-4),f(1E+4)])))
                nu_new = fzero(f, [1E-4,1E+4]);
            else
                factor = 2;
                Mm = [nu/factor, nu*factor];
                while f(Mm(1))*f(Mm(2)) >= 0
                    Mm(1) = Mm(1) / factor;
                    Mm(2) = Mm(2) * factor;
                end
                if(and(prod(isfinite(Mm)), prod(isfinite(f(Mm)))))
                    nu_new = fzero(f, Mm);
                else
                    break
                end
            end
        else
            factor = 2;
            Mm = [nu/factor, nu*factor];
            while f(Mm(1))*f(Mm(2)) >= 0
                Mm(1) = Mm(1) / factor;
                Mm(2) = Mm(2) * factor;
            end
            if(and(prod(isfinite(Mm)), prod(isfinite(f(Mm)))))
                nu_new = fzero(f, Mm);
            else
                break
            end
        end
        nu = nu_new;
        mu = mu_new;
        lambda = lambda_new;
    end
else
    mu = mu_init;
    lambda = lambda_init;
    nu = nu_init;
    
    %   Robust mixture modelling using the t distribution
    %   D. Peel, G. J. McLachlan
    %   Statistics and Computing (2000)
    
    u = (nu+1)*ones(size(data))./(nu + lambda*(data-mu).^2);
    mu_new = sum(W.*u.*data)/(u'*W);
    lambda_new = (sum(W.*u.*(data - mu_new).^2)/sum(W))^-1;
%     f = @(nu) log(nu) - psi(nu/2) +1 - log(nu+1) + psi((nu+1)/2) + (W'*(log(u)-u))/sum(W);
    
    opts = optimoptions(@fmincon,'Display','off');
    opts.Algorithm = 'sqp';
    nu_new = fmincon(@(nu)0,nu,[],[],[],[],0,[],@(nu) computeNewNu(nu, W, u),opts);
%     if(prod(isfinite([f(1E-4),f(1E+4)])))
%         nu_new = fzero(f, [1E-4,1E+4]);
%     else
%         factor = 2;
%         Mm = [nu/factor, nu*factor];
%         while f(Mm(1))*f(Mm(2)) >= 0
%             Mm(1) = Mm(1) / factor;
%             Mm(2) = Mm(2) * factor;
%         end
%         if(and(prod(isfinite(Mm)), prod(isfinite(f(Mm)))))
%             nu_new = fzero(f, Mm);
%         else
%             nu_new = nu;
% %             error('Impossible to estimate nu of the Student-t distribution');
%         end
%     end
    
    mu = mu_new; lambda = lambda_new; nu = nu_new;
end
end

function [c,ceq] = computeNewNu(nu, W ,u)

c = []; % no nonlinear inequality
ceq = log(nu) - psi(nu/2) +1 - log(nu+1) + psi((nu+1)/2) + (W'*(log(u)-u))/sum(W); % the fsolve objective is fmincon constraints
end