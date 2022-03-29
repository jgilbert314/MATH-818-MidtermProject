function [ D_op ] = calcDispOp(alpha, a, a_dag)

% Input:
%   alpha - complex number
%   a - matrix
%   a_dag - matrix (complex transpose of a)


X = alpha*a - conj(alpha)*a_dag;
D_op = expm(X);

end