function [ alpha_vec ] = calcDispState(alpha, n_vec)

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