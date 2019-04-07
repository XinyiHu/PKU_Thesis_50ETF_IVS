function [prices, alphas] = hestonprice(type,S,K,T,r,q,v0,theta,rho,kappa,sigma,alphas)

% Function to calculate European option price based on Heston Model
% Computation method: Lord, Kahl (1993) Optimal Fourier inversion in 
% semi-analytical option pricing

% Input: (PC till q can be vectorized)
% type: 'EuropeanCall' or 'EuropeanPut'
% S: Spot
% K: Strike
% T: Year to maturity
% r: interest rate
% q: dividend rate
% v0: initial variance
% theta: long run mean variance
% kappa: mean reversion speed of  volatility
% sigma: volatility of volatility
% rho: correlation between return volatility
% alpha: alpha can be a vector supplied by the user, otherwise the
%       function attempts to find a payoff-dependent optimal alpha


prices = NaN(size(type,1),1);
F = S .* exp((r-q).*T);

if(~exist('alphas','var'))
    alphas = NaN(size(type,1),1);
elseif(numel(alphas)==1)
    alphas = repmat(alphas,numel(S),1);
end

alpha0=0.75;

for i=1:size(type,1)
    if (isnan(alphas(i)))
        try
            alphas(i) = fzero(@(a) psi(a,K(i),F(i),kappa,theta,rho,sigma,T(i),v0), alpha0);
        catch
           alphas(i) = alpha0;
        end
    end
    prices(i) = Ralpha(F(i),K(i),alphas(i)) + 1/pi*integral(@(x) phi(x,K(i),alphas(i),F(i),kappa,theta,rho,sigma,T(i),v0) , 0, Inf);
    if strcmp(type(i),'EuropeanPut')
        prices(i) = prices(i) + K(i)*exp(-r(i)*T(i)) - S(i)*exp(-q(i)*T(i));
    end
end

    
end


function p = psi(alpha,K,F,kappa,theta,rho,sigma,T,v0)
    k = log(K);
    p = -alpha*k + 0.5*log(phi(-(alpha+1)*1i,K,alpha,F,kappa,theta,rho,sigma,T,v0)^2);
end

function r = Ralpha(F,K,alpha)
    r = F*(alpha<=0) - K*(alpha<=-1) - 0.5*(F*(alpha==0)-K*(alpha==-1));
end

function y = phi(v,K,alpha,F,kappa,theta,rho,sigma,T,v0)
    k = log(K);
    y = real(exp(-1i*(v-1i*alpha)*k) .* (cf(v-1i*(alpha+1),F,kappa,theta,rho,sigma,T,v0) ./ (-(v-1i*(alpha+1)).*(v-1i*alpha))));
end

function c = cf(u,F,kappa,theta,rho,sigma,T,v0)
    f = log(F);
    c = exp(1i*u*f+ A(u,kappa,theta,rho,sigma,T) + Bv(u,rho,sigma,kappa,T)*v0);
end

function b = Bv(u,rho,sigma,kappa,T)
    b = ((beta(u,rho,sigma,kappa)-D(u, rho, sigma, kappa)) .* (1-exp(-D(u, rho, sigma, kappa)*T))) ./ (sigma.^2*(1-G(u, rho, sigma, kappa).*exp(-D(u, rho, sigma, kappa)*T)));
end

function a = A(u,kappa,theta,rho,sigma,T)
    a = (kappa*theta*((beta(u,rho,sigma,kappa)-D(u,rho,sigma,kappa))*T - 2*log(phi2(u,rho,sigma,kappa,T))))/sigma.^2;
end

function p = phi2(u,rho,sigma,kappa,T)
    p = (G(u,rho,sigma,kappa) .* exp(-D(u, rho, sigma, kappa)*T)-1) ./ (G(u, rho, sigma, kappa)-1);
end

function g = G(u,rho,sigma,kappa)
    g = (beta(u,rho,sigma,kappa)-D(u, rho, sigma, kappa)) ./ (beta(u,rho,sigma,kappa)+D(u, rho, sigma, kappa));
end

function d = D(u,rho,sigma,kappa)
    d = sqrt(beta(u,rho,sigma,kappa).^2 - 4*alphahat(u)*gamma(sigma));
end

function a = alphahat(u)
    a = -0.5*u.*(1i+u);
end

function b = beta(u,rho,sigma,kappa)
    b = kappa-rho*sigma*u*1i;
end

function y = gamma(sigma)
    y = 0.5*sigma.^2;
end
