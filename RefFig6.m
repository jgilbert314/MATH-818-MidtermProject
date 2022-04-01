
N_n = 1e2 + 1;
N_g = 1;
n_g = linspace(-N_g, N_g, N_n);
% n_g = zeros(1, N_n);


test = zeros(N, N_n);
cos_mat = ( expm(1i*vPhi_hat) + expm(-1i*vPhi_hat) )/2;
for itr = 1:N_n
    H_hat = 4*E_C*(n_hat - n_g(itr))^2 - E_J*cos_mat;
    
    these_eig = sort( eig(H_hat) );
    test(:, itr) = these_eig - these_eig(1);
end
test = test/hbar/2/pi * 1e-9;

plot(w_r, abs(test(:, :))', '.');


%% 

[V, D] = eig(a_dag);