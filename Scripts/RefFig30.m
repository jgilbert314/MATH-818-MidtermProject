%% Calculates eigenvalues of system
% Meant to reacreate Fig 30

% N defined externally
if (N < 5) % Figure uses at least 5 energy levels
   error('Requires N > 5'); 
end

% Ladder operators
[a, a_dag] = calcOscLadderOp(N);


% Oscillator Basis
[ Phi_hat, Q_hat ] = calcOscBasis(a, a_dag, Z_r, hbar);

% Test of scaling. Should evaluate to identity
% [Phi_hat, Q_hat] = i*hbar
testOsc = (Phi_hat*Q_hat - Q_hat*Phi_hat)/1i/hbar;


% Transmon Basis (depends on oscillator basis) 
[ n_hat, vPhi_hat ] = calcTransBasis(Phi_hat, Q_hat, hbar, e);
% Ladder operators
[b, b_dag] = calcTransLadderOp(a, a_dag, Z_r, E_C, E_J, hbar, e);

% Test of scaling. Should evaluate to identity
testTrans = (vPhi_hat*n_hat - n_hat*vPhi_hat)/1i;


% Parameters defined in Figure caption
E_C1 = 2*pi*hbar * 317e6;
E_C2 = 2*pi*hbar * 297e6;
E_J1 = 2*pi*hbar * 28.48e9;
E_J2 = 2*pi*hbar * 42.34e9;
g1 = 2*pi * 199e6;
g2 = 2*pi * 190e6;

% Bounds of parametric sweep
w_q1 = linspace(1, 50, 500)*1e9 * 2*pi;
E_J1 = (hbar*w_q1 + E_C1).^2 / 8/E_C1;


% Run Parametric sweep
param_vec = E_J1;
hamilHandle = ...
    @(E_J1) calcOscCoupledTransHamiltonian(g1, g2, E_C1, E_J1, E_C2, E_J2, a, a_dag, w_r, Z_r, hbar, e);

[ eig_arr ] = sweepHamilEig(param_vec, hamilHandle);

% Plot Results
tran_freq = w_q1/2/pi * 1e-9; % Convert to GHz
eig_freq = eig_arr/hbar/(2*pi)^2 * 1e-9; % Convert to GHz
plot(tran_freq, eig_freq(1:6, :)');
xlabel('Second Qubit Frequency [GHz]');
ylabel('Eigenfrequency');
% xlim([10, 20]);
ylim([0, 20]);
grid on;
title(['N = ', num2str(N), '    Z_r = ', num2str(Z_r)]);
legend('2nd Qubit', 'Resonator', '1st Qubit');
drawnow;