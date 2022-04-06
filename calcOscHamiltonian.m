function [ H_s ] = calcOscHamiltonian(a, a_dag, w_r, hbar)
% Summary:
%   Defines the Hamiltonian for a harmonic oscillator
%       H_s = hbar*w_r*(a_dag*a - 1/2)
% 
% Input:
%   a, a_dag - creation and ahnihilation operators (square matrices)
%               a_dag = ctranspose(a);
%   w_r - resonant frequency of microwave cavity
%               From paper (p.5): w_r ~ 2*pi * 8 GHz
%   kappa - damping coefficient
%   hbar - reduced Planck's constant
%           hbar ~ 1.0545718 × 10-34 m2 kg / s
% 
% Output:
%   H_s - square matrix representing the Hamiltionian

% Author: Jason Gilbert
% Date: March 15, 2022
% Version: V00
% Last Updated: N/A

H_s = hbar*w_r*(a_dag*a - 1/2);

end