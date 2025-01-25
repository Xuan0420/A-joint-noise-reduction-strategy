function Entropy = BaoluoEntropy(y)
% y - IMF components after VMD decomposition
% Hilbert transform to obtain the analytic signal
yh = hilbert(y);    % MATLAB function returns a complex signal, the analytic signal
% Imaginary part
xi = imag(yh); 
% Real part
xr = real(yh);
% Envelope signal
A = sqrt(xr.^2 + xi.^2);
% Probability distribution sequence
p = A ./ (sum(A));
% Envelope entropy
Entropy = -sum(p .* log(p));  
