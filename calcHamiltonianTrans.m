function [ H ] = calcHamiltonianTrans(n_hat, Vphi_hat, E_C, E_J)

n_sq = n_hat*n_hat;
vphi_sq = Vphi_hat*Vphi_hat;

H = 4*E_C*n_sq + E_J*vphi_sq*(1/2 - vphi_sq/32);

end