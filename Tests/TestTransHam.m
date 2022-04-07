%% Initialization
clear;
close all;
DefineConstants; % Run initialization script


%%
N = 5;

% Calc ladder operators
[a, a_dag] = calcOscLadderOp(N);
[b, b_dag] = calcTransLadderOp(a, a_dag, Z_r, E_C, E_J, hbar, e);


% Oscillator Basis
[Phi_hat, Q_hat] = calcOscBasis(a, a_dag, Z_r, hbar);
[n_hat, vPhi_hat] = calcTransBasis(Phi_hat, Q_hat, hbar, e);

% Test of scaling. Should evaluate to identity
% [Phi_hat, Q_hat] = i*hbar
testOsc = (Phi_hat*Q_hat - Q_hat*Phi_hat)/1i/hbar;
testTrans = (n_hat*vPhi_hat - vPhi_hat*n_hat)/1i;


% Density matrix
[psi_base, ~] = eig(n_hat);
psi = psi_base(:, 1);
rho = kron(psi', psi); % Density matrix


%% Initialize Input Structure

InputStruct.a = a;
InputStruct.a_dag = a_dag;
InputStruct.kappa = kappa;
InputStruct.n_k = n_k;

disp('Running Single Transmon');
InputStruct.hamilHandle = ...
    @ (t) calcTransHamiltonian(n_hat, vPhi_hat, E_C, E_J);



%% Calculate Solution

rho_0 = reshape(rho, numel(rho), 1);
% funcHandle = @(t, y) masterEq(y, a, a_dag, w_r, kappa, n_k, hbar);
funcHandle = @(t, y) masterEq(t, y, InputStruct);
[t_out, y_out] = ode45(funcHandle, [0, 10], rho_0);



%% Plot Solution

diag_inds = 1 + [0:N]*(N+2);
y_diag = y_out(:, diag_inds); % Diagonal elements of rho

figure(1);
plot(t_out, abs(y_diag))
xlabel('Time');
ylabel('Magnitude |\rho|');
title('Magnitude of Diagonal Elements');
grid on;

figure(2);
imagesc(1:numel(rho), t_out, abs(y_out));
ylabel('Time');
xlabel('Matrix Index');
title('Density Matrix Magnitudes');
colorbar;


%% Check trace

N_t = length(t_out);
rho_tr = zeros(1, N_t);
for itrT = 1:N_t
   this_rho = reshape(y_out(itrT, :), [N+1, N+1]);
   rho_tr(itrT) = trace( (this_rho)^2 );
end

plot(t_out, rho_tr);