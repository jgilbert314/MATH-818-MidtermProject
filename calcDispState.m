function [ alpha_vec ] = calcDispState(alpha, n_vec)
% WARNING: does not work yet, do not use
% 
% Calculates coherent state defined by Eq.87
% 
% Input:
%   alpha - phase parameter (complex value)
%   n_vec - not sure yet
% 
% Output:
%   alpha_vec - vector of coherent state

% Author: Jason Gilbert
% Date: ??
% Version: N/A
% Last Updated: N/A

tol = 1e-12;
alpha_vec = zeros(length(n_vec), 1);
ind = 0;

while (ind < 100)
    term = alpha^ind/sqrt(factorial(ind))*n_vec;
    alpha_vec = alpha_vec + term;
    
    ind = ind+1;
end

alpha_vec = exp(-abs(alpha)^2/2)*alpha_vec;

end