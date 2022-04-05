%% Calculates eigenvalues of system
% Meant to reacreate Fig 30

E_C1 = 2*pi*hbar * 317e6;
E_C2 = 2*pi*hbar * 297e6;
E_J1 = 2*pi*hbar * 28.48e9;
E_J2 = 2*pi*hbar * 42.34e9;
g1 = 2*pi * 199e6;
g2 = 2*pi * 190e6;

w_q2 = linspace(1, 50, 500)*1e9 * 2*pi;
E_J2 = (hbar*w_q2 + E_C2).^2 / 8/E_C2;


param_vec = E_J2;

hamilHandle = ...
    @(E_J2) calcOscCoupledTransmons(g1, g2, E_C1, E_J1, E_C2, E_J2, a, a_dag, w_r, Z_r, hbar, e);

[ eig_arr ] = sweepHamilEig(param_vec, hamilHandle);
plot(param_vec, eig_arr');


tran_freq = w_q2/2/pi * 1e-9;
eig_freq = eig_arr/hbar/(2*pi)^2 * 1e-9;
plot(tran_freq, eig_freq');
xlabel('Second Qubit Frequency [GHz]');
ylabel('Eigenfrequency');
% xlim([10, 20]);
ylim([0, 20]);
grid on;

