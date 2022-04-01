function [ H ] = calcHamiltonianTransApprox(b, b_dag, E_C, E_J)
% An approximation of the transmon Hamiltonian, as defined in Eq.(27), p.10
% 
% Input:
%   b, b_dag - transmon ladder operators
%   E_C, E_J - Charging energy, Josephson energy
% 
% Output:
%   H - hamiltonian (square matrix)

% Author: Jason Gilbert
% Date: March 31, 2022
% Version:
% Last Updated: N/A

coeff = sqrt(8*E_C*E_J) - E_C;

% H = coeff*b_dag*b - (E_C*b_dag^2*b^2)/2;
H = coeff*b_dag*b - (E_C*b_dag*b_dag*b*b)/2;

end