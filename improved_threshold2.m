function [w_hat] = improved_threshold2(w, lambda, alpha, beta)
% IMPROVED_THRESHOLD Calculates the threshold based on the improved soft-threshold function
%   w: Input data vector
%   lambda: Threshold parameter
%   alpha: Exponent parameter
%   beta: Exponent parameter

% Calculate the absolute value and sign
abs_w = abs(w);
sgn_w = sign(w);

% Calculate the threshold based on the improved threshold function
w_hat = zeros(size(w));
idx = (abs_w >= lambda);
w_hat(idx) = sgn_w(idx) .* (abs_w(idx)-(abs(lambda) ./ ((abs_w(idx).^beta - abs(lambda).^beta + 1).^(1/beta) .* exp((abs_w(idx) - abs(lambda)).^(1/alpha)))));

end
