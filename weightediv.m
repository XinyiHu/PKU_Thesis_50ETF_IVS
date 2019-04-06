function wiv = weightediv(iv, volume)

% Function to calculate volume weighted implied volatility of multiple options
% Input:
% iv: N*1 array of implied volatility
% volume: N*1 array of option volume

wiv = sum(iv .* volume) / sum(volume);

end