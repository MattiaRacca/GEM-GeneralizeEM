function [ y ] = Student( mu, lambda, nu, x )
%STUDENT Student-t probability distribution
%   http://www.inference.phy.cam.ac.uk/cs482/publications/scheffler2008derivation.pdf

y = gamma((nu+1)/2)/gamma(nu/2)*sqrt((lambda)/(pi*nu))*(1 + (lambda*(x-mu).^2)/(nu)).^(-(nu+1)/2);

end

