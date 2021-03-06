function [ H ] = calcOscCoupledTransHamiltonian(g1, g2, E_C1, E_J1, E_C2, E_J2, a, a_dag, w_r, Z_r, hbar, e)
% Calculates the Hamiltonian of two transmons coupled via a harmonic
% oscillator.
%   From Eq.(138) of paper (p.42)
% 
% Input:
%   g1, g2 - coupling constants of each transmon (should be proportional to
%            w_q or w_r)
%   E_C, E_J - charging, Josephson energies for each transmon
%   a, a_dag - ladder operators for Harmonic oscillator
%   omega_r - resonant frequency of oscillator
%   Z_r - characteristic impedance of oscillaotr
%   hbar, e - physical constants
% 
% Output:
%   H - square matrix representing Hamiltonian

% Author: Jason Gilbert
% Date: March 31, 2022
% Version: N/A
% Last Updated: April 04, 2022
%   Updated to calculate transmon ladder operators instead of taking them
%   as inputs. This ensures that the ladder operators correspond to the
%   parameters of the transmon

[b1, b1_dag] = calcTransLadderOp(a, a_dag, Z_r, E_C1, E_J1, hbar, e);
[b2, b2_dag] = calcTransLadderOp(a, a_dag, Z_r, E_C2, E_J2, hbar, e);

% Hq = hbar*w_q*b_dag*b - E_C*b_dag^2*b^2/2
%   Eq.(27)
%   hbar*w_q = sqrt(8*E_C*E_J) - E_C
H1 = calcTransApproxHamiltonian(b1, b1_dag, E_C1, E_J1);
H2 = calcTransApproxHamiltonian(b2, b2_dag, E_C2, E_J2);

% Oscillator Term:
%   hbar*w_r*a_dag*a
T0 = hbar*w_r*a_dag*a;

% Interaction Terms
%   hbar*g*(a_dag*b + a*b_dag)
T1 = hbar*g1*(a_dag*b1 + a*b1_dag);
T2 = hbar*g2*(a_dag*b2 + a*b2_dag);

% Full Hamiltonian
H = H1 + H2 + T0 + T1 + T2;

end