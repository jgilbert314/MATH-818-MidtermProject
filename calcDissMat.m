function [ D_mat ] = calcDissMat(A, B)
% Calculates the dissipative terms in the master equation (Eq 70)
%   D_mat = D[A]B = A*B*A_dag - 1/2*(C*B + B*C)
%   C = A_dag*A

% Author: Jason Gilbert
% Date: March 10, 2022
% Version: V00
% Last Updated: N/A

A_dag = ctranspose(A); % Hermitian transpose
C = A_dag*A;           % Intermediate operator

% Term 1
term_1 = A*B*A_dag;

% Term 2 (anti-commutator)
term_2 = C*B + B*C;

D_mat = term_1 - term_2/2;


end