function [ H ] = calcHamiltonianDrive(a, a_dag, w_r, kappa, w_d, phi_d, ampHandle, t, hbar)
% Summary:
%   Defines the Hamiltonian for a harmonic oscillator with a driving force
% 
% Input:
%   a, a_dag - creation and ahnihilation operators (square matrices)
%               a_dag = ctranspose(a);
%   w_r - resonant frequency of microwave cavity
%               From paper (p.5): w_r ~ 2*pi * 8 GHz
%   kappa - damping coefficient
%   w_d - drive frequency
%   phi_d - drive phase
%   ampHandle - function handle defining the drive amplitude A(t)
%   hbar - reduced Planck's constant
%           hbar ~ 1.0545718 × 10e-34 m2 kg / s
% 
% Output:
%   H - square matrix representing the Hamiltionian

% Author: Jason Gilbert
% Date: March 16, 2022
% Version: V00
% Last Updated: N/A

H_s = hbar*w_r*(a_dag*a - 1/2);

eps = 1i*sqrt(kappa)*ampHandle(t);
wave_term = eps*exp( -1i*(w_d*t - phi_d) );
H_d = hbar*( a*wave_term + a_dag*conj(wave_term) );

H = H_s + H_d;

end