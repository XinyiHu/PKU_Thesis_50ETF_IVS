function sigma = bsmiv(type,price,S,K,T,r,q)

% Function to calculate the implied volatility of the BSM model
% Approximation: Brenner and Subrahmanyam (1988) 
% Solution finding: fzero
% Input:
% price:option price
% S: underlying spot price
% K: strike price
% T: year to maturity
% t: date of spot price
% r: interest rates
% q: dividend rates


fun = @(price_,type_,S_,K_,T_,r_,q_,sigma_) price_ - bsmprice(type_,S_,K_,T_,r_,q_,sigma_);

if (size(K,1) ~= size(T,1) || size(K,1) ~= size(type,1) || size(K,1) ~= size(S,1) || size(price,1) ~= size(S,1))
   disp ('type,price,S,K,T,r,q must have same dimensions');
   return
end

sigma = NaN(size(type,1),1);

for i=1:size(K,1)
    sigma_star = sqrt(2*pi/T(i)) * (price(i)/S(i));
    try
        sigma(i) = fzero(@(x)fun(price(i),type(i),S(i),K(i),T(i),r(i),q(i),x), sigma_star);
    catch
        sigma(i) = fzero(@(x)fun(price(i),type(i),S(i),K(i),T(i),r(i),q(i),x), 0.3);
    end
end

end

