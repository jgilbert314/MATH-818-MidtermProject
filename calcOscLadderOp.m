function [ a, a_dag ] = calcOscLadderOp(N)
% Calculates ladder operators for a harmonic oscillator
% 
% Input:
%   N - number of allowed states
% 
% Output:
%   a, a_dag - ladder operators (N+1) by (N+1) matrices

% Author: Jason Gilbert
% Date: April 06, 2022
% Version: V00
% Last Updated: N/A

a_base = sqrt( 1:N ); 
a = diag(a_base, 1);
a_dag = ctranspose(a);

end