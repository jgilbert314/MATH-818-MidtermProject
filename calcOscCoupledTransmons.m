function [ H ] = calcOscCoupledTransmons(b1, b1_dag, b2, b2_dag, g1, g2, E_C1, E_J1, E_C2, E_J2, a, a_dag, w_r, hbar)
% From Eq.(138) of paper (p.138)
% 
% Input:
%   a, a_dag
%   b1, b1_dag
%   b2, b2_dag
%   omega_r
%   g1, g2 - should be proportional to w_q or w_r

% Hq = hbar*w_q*b_dag*b - E_C*b_dag^2*b^2/2
%   Eq.(27)
%   hbar*w_q = sqrt(8*E_C*E_J) - E_C
H1 = calcHamiltonianTransApprox(b1, b1_dag, E_C1, E_J1);
H2 = calcHamiltonianTransApprox(b2, b2_dag, E_C2, E_J2);

% Oscillator Term:
%   hbar*w_r*a_dag*a
T0 = hbar*w_r*a_dag*a;

% Interaction Terms
%   hbar*g*(a_dag*b + a*b_dag)
T1 = hbar*g1*(a_dag*b1 + a*b1_dag);
T2 = hbar*g2*(a_dag*b2 + a*b2_dag);

H = H1 + H2 + T0 + T1 + T2;

end