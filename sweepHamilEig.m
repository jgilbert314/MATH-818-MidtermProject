function [ eig_arr ] = sweepHamilEig(param_vec, hamilHandle)

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