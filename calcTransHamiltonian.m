function [ H ] = calcTransHamiltonian(n_hat, vPhi_hat, E_C, E_J)
% Approximate Hamiltonian for transmon, as define in Eq.(24)
% 
% Input:
%   n_hat - square matrix
%   vPhi_hat - square matrix
%   E_C - constant
%   E_J - constant
% 
% Output:
%   H - square matrix

% Author: Jason Gilbert
% Date: March ??, 2022

% TODO:
% Test full transmon Hamiltonian:
% cos_mat = ( expm(1i*vPhi_hat) + expm(-1i*vPhi_hat) )/2;
% H_hat = 4*E_C*(n_hat - n_g).^2 - E_J*cos_mat;

n_sq = n_hat^2;
vphi_sq = vPhi_hat^2;

H = 4*E_C*n_sq + E_J*vphi_sq*(1/2 - vphi_sq/32);

end