function [ D_op ] = calcDispOp(alpha, a, a_dag)
% Calculates displacement operator:
%   D_op = exp(X) 
%       where X = alpha*a - conj(alpha)*a_dag
%   Matrix exponentiation calculated via D_op = expm(X)
% 
% Input:
%   alpha - complex number
%   a - matrix
%   a_dag - matrix (complex transpose of a)
% 
% Output:
%   D_op - matrix of same size as a and a_dag

% Author: Jason Gilbert
% Date: March 29, 2022

X = alpha*a - conj(alpha)*a_dag;
D_op = expm(X);

end