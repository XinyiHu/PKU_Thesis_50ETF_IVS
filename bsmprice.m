function price = bsmprice(type,S,K,T,r,q,sigma)

% Function to calculate Black-Scholes-Merton Price for european options
% Input:
% type: 'EuropeanCall' or 'EuropeanPut'
% S: underlying spot price
% K: strike price
% T: year to maturity
% r: interest rate
% q: dividend rate
% sigma: volatility of returns of underlying

    
d1 = (log(S./K) + ((r-q) + 0.5*sigma.^2) .* T) ./ (sigma .* sqrt(T));
d2 = d1 - sigma .* sqrt(T);

price = NaN(size(type,1),1);

for i=1:length(type)
    if strcmp(type(i), 'EuropeanCall')
        price(i) = S(i)*normcdf(d1(i)) - K(i)*exp(-(r(i)-q(i))*T(i))*normcdf(d2(i));
    elseif strcmp(type(i), 'EuropeanPut')
        price(i) = normcdf(-d2(i))*K(i)*exp(-(r(i)-q(i))*T(i)) - normcdf(-d1(i))*S(i);
    else
        disp('Invalid option type'); 
    end
end

end

