function [b, b_dag] = calcTransLadderOp(a, a_dag, Z_r, E_C, E_J, hbar, e)
% Calculates ladder operators for the transmon (Josephson junction)
% Formula obtained from Eq.(25-26) in reference
% 
% Inputs:
%   a, a_dag - harmonic ladder operators
%   Z_r - impedance
%   E_C, E_J - charging, Josephson energies
%   hbar, e - physical constants
% 
% Output:
%   b, b_dag - ladder operators

% Author: Jason Gilbert
% Date: March 31, 2022
% Version: 
% Last Updated: April 06, 2022
%   Bases calculated via functions

% Oscillator Basis
[Phi_hat, Q_hat] = calcOscBasis(a, a_dag, Z_r, hbar);

% Transmon Basis
vPhi_0 = (2*E_C/E_J)^(1/4);
nHat_0 = 1i/2*(E_J/2/E_C)^(1/4); 
[n_hat, vPhi_hat] = calcTransBasis(Phi_hat, Q_hat, hbar, e);


% Ladder operators
b = (vPhi_hat/vPhi_0 - n_hat/nHat_0)/2;
b_dag = (vPhi_hat/vPhi_0 + n_hat/nHat_0)/2;


end