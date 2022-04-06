function [ n_hat, vPhi_hat ] = calcTransBasis(Phi_hat, Q_hat, hbar, e)
% Calculates operators corresponding to canonical coordinates of transmon
%   [ Phi_hat, Q_hat ] = calcOscBasis(a, a_dag, Z_r, hbar);
%   n_hat = Q_hat/2/e
%   vPhi_hat = 2*e/hbar*Phi_hat
% 
% Input:
%   Phi_hat - square matrix
%   Q_hat - square matrix
%   hbar - reduced Planck's constant
%   e - charge of electron
% 
% Output:
%   n_hat - square matrix
%   vPhi_hat - square matrix

% Author: Jason Gilbert
% Date: ??
% Version: N/A
% Last Updated: N/A

% Transmon Basis (depends on oscillator basis) 
n_hat = Q_hat/2/e; % Defined below Eq.(22)
vPhi_hat = 2*e/hbar*Phi_hat; % Phi_0 defined below Eq.(19)

end