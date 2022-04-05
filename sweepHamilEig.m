function [ eig_arr ] = sweepHamilEig(param_vec, hamilHandle)
% Calculates the eigenvalues of a Hamilonian as a function of an arbitrary
% parameter.
% 
% Input:
%   param_vec - vector of parameter values
%   hamilHandle - a function handle specifying the Hamiltonian
% 
% Output:
%   eig_arr - an NxM matrix, where N is the number of eigenvalues, and M is
%             the number of parameters. Sorted in ascending order with N.

% Author: Jason Gilbert
% Date: March 31, 2022
% Version: N/A
% Last Updated: N/A


N_p = length(param_vec);
H_0 = hamilHandle(param_vec(1));

eig_arr = zeros(length(H_0), N_p);

eig_arr(:, 1) = sort( eig(H_0) );
for itr = 2:N_p
    this_H = hamilHandle(param_vec(itr));
    these_eig = eig(this_H);
    
    eig_arr(:, itr) = sort(these_eig);
end

end