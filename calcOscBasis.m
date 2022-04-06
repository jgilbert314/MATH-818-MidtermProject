function [ Phi_hat, Q_hat ] = calcOscBasis(a, a_dag, Z_r, hbar)
% Calculates operators corresponding to canonical coordinates of oscillator
% 	Phi_0 = sqrt(hbar*Z_r/2);
% 	Q_0 = sqrt(hbar/2/Z_r);
% 
%	Phi_hat = Phi_0*(a_dag + a);
%	Q_hat = 1i*Q_0*(a_dag - a);
% 
% Input: 
%   a, a_dag - ladder operators
%   Z_r - characteristic impedance
%   hbar - reduced Planck's constant
% 
% Output:
%   Phi_hat - flux operator (square matrix)
%   Q_hat - charge operator (square matrix)

% Author: Jason Gilbert
% Date: ??
% Version: N/A
% Last Updated: N/A

Phi_0 = sqrt(hbar*Z_r/2);
Q_0 = sqrt(hbar/2/Z_r);

Phi_hat = Phi_0*(a_dag + a);
Q_hat = 1i*Q_0*(a_dag - a);

end