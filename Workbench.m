%% Initializaiton
clear;

psi = [1 ; 0 ; 0 ; 0]; % State vectors
a_base = sqrt( 1:length(psi)-1 ); 
a = diag(a_base, 1);
a_dag = ctranspose(a);

rho = kron(psi', psi); % Density matrix

% Constants
hbar = 1.0545718e-34; % Reduced Planck's constant [m^2*kg/s]
n_k = 1.8;            % Constant from Bose-Einstein distribution
kappa = 0;            % Photon decay rate
w_r = 2*pi*8e9;       % Cavity frequency

Z_r = 50;             % Transmission Line Impedance (Ohms)


% Basis
Phi_0 = sqrt(hbar*Z_r/2);
Q_0 = sqrt(hbar/2/Z_r);

Phi_hat = Phi_0*(a_dag + a);
Q_hat = 1i*Q_0*(a_dag - a);

% Test of scaling. Should evaluate to identity
% [Phi_hat, Q_hat] = i*hbar
test = (Phi_hat*Q_hat - Q_hat*Phi_hat)/1i/hbar;


%% Initialize Input Structure

w_d = w_r;
phi_d = pi*1.15;
ampHandle = @(t) 10/hbar;

InputStruct.a = a;
InputStruct.a_dag = a_dag;
InputStruct.kappa = kappa;
InputStruct.n_k = n_k;
InputStruct.hamilHandle = ...
    @ (t) calcHamiltonianDrive(a, a_dag, w_r, kappa, w_d, phi_d, ampHandle, t, hbar);

%     @ (t) calcHamiltonianOsc(a, a_dag, w_r, hbar)
%     @ (t) calcHamiltonianDrive(a, a_dag, w_r, kappa, w_d, phi_d, ampHandle, t, hbar);
%     @ (t) calcHamiltonianTrans(Phi_hat, Q_hat, 1, 1);

%% Calculate Solution


rho_0 = reshape(rho, numel(rho), 1);
% funcHandle = @(t, y) masterEq(y, a, a_dag, w_r, kappa, n_k, hbar);
funcHandle = @(t, y) masterEq(t, y, InputStruct);
[t_out, y_out] = ode45(funcHandle, [0, 1], rho_0);



%% Plot Solution

y_diag = y_out(:, :); % Diagonal elements of rho

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