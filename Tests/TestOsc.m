%% Initializaiton
clear;
close all;
DefineConstants; % Define hbar and Z_r

N = 5;
psi = [0:N]; % State vectors


% Calc ladder operators
[a, a_dag] = calcOscLadderOp(N);

% Calc basis operators
[Phi_hat, Q_hat] = calcOscBasis(a, a_dag, Z_r, hbar);

% Test of scaling. Should evaluate to the following identities
% [Phi_hat, Q_hat] = i*hbar
testOsc = (Phi_hat*Q_hat - Q_hat*Phi_hat)/1i/hbar;


%% Calculate Errors

% Self-adjointness
err_adj = sum(sum( (a - ctranspose(a_dag)).^2 ));

% Commutation relation
I_mat = eye(N-1);
osc_mat = testOsc(1:N-1, 1:N-1);

err_osc = sum(sum( (osc_mat - I_mat).^2 ));

