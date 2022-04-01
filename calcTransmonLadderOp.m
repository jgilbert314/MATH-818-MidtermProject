function [b, b_dag] = calcTransmonLadderOp(a, a_dag, Z_r, E_C, E_J, hbar, e)
% 
% 
% Inputs:
%   Z_r - impedance
%   E_C, E_J - energies
%   a, a_dag - harmonic ladder operators
%   hbar, e - physical constants
% 
% Output:
%   b, b_dag - ladder operators

% Author: Jason Gilbert
% Date: March 31, 2022
% Version: 
% Last Updated: N/A

% Oscillator Basis
Phi_0 = sqrt(hbar*Z_r/2);
Q_0 = sqrt(hbar/2/Z_r);

Phi_hat = Phi_0*(a_dag + a);
Q_hat = 1i*Q_0*(a_dag - a);


% Transmon bais
vPhi_0 = (2*E_C/E_J)^(1/4);
nHat_0 = 1i/2/vPhi_0; 

n_hat = Q_hat/2/e;
vPhi_hat = 2*e/hbar*Phi_hat;


% Ladder operators
b = (vPhi_hat/vPhi_0 - n_hat/nHat_0)/2;
b_dag = (vPhi_hat/vPhi_0 + n_hat/nHat_0)/2;


end