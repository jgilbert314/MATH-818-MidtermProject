%% Initializaiton
clear;
close all;
DefineConstants; % Run initialization script

N = 5;
psi = [0:N]; % State vectors


% Calc ladder operators
[a, a_dag] = calcOscLadderOp(N);
[b, b_dag] = calcTransLadderOp(a, a_dag, Z_r, E_C, E_J, hbar, e);

% Calc basis operators
[Phi_hat, Q_hat] = calcOscBasis(a, a_dag, Z_r, hbar);
[n_hat, vPhi_hat] = calcTransBasis(Phi_hat, Q_hat, hbar, e);

% Test of scaling. Should evaluate to the following identities
% [Phi_hat, Q_hat] = i*hbar
testOsc = (Phi_hat*Q_hat - Q_hat*Phi_hat)/1i/hbar;
% [n_hat, vPhi_hat] = -1
%   Zagoskin, p.68
testTrans = (n_hat*vPhi_hat - vPhi_hat*n_hat)/-1;
% [vPhi_hat, Q_hat] = 2i*e
%   Martinis, p.6
testTrans2 = (vPhi_hat*Q_hat - Q_hat*vPhi_hat)/2i/e;


% Test self-adjointness
err_a = sum(sum( (a - ctranspose(a_dag)).^2 ));
err_b = sum(sum( (b - ctranspose(b_dag)).^2 ));