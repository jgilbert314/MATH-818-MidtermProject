%% Initializaiton
clear;
DefineConstants; % Run initialization script

psi = [0:5]; % State vectors
N = length(psi);

a_base = sqrt( 1:length(psi)-1 ); 
a = diag(a_base, 1);
a_dag = ctranspose(a);


% Oscillator Basis
Phi_0 = sqrt(hbar*Z_r/2);
Q_0 = sqrt(hbar/2/Z_r);

Phi_hat = Phi_0*(a_dag + a);
Q_hat = 1i*Q_0*(a_dag - a);

% Test of scaling. Should evaluate to identity
% [Phi_hat, Q_hat] = i*hbar
test = (Phi_hat*Q_hat - Q_hat*Phi_hat)/1i/hbar;


% Transmon Basis (depends on oscillator basis) 
n_hat = Q_hat/2/e;
vPhi_hat = 2*pi/Phi_0*Phi_hat;

% Ladder operators
vPhi_0 = (2*E_C/E_J)^(1/4);
nHat_0 = 1i/2/vPhi_0; 
b = (vPhi_hat/vPhi_0 - n_hat/nHat_0)/2;
b_dag = (vPhi_hat/vPhi_0 + n_hat/nHat_0)/2;


% Density matrix
psi_0 = zeros(1, N);
psi_1 = psi_0;
psi_1(1) = 1;
psi_2 = psi_0;
psi_2(2) = 1;
rho = kron(psi', psi); % Density matrix


%% Initialize Input Structure

w_d = w_r;
phi_d = pi*1.15;
ampHandle = @(t) 10/hbar;

InputStruct.a = a;
InputStruct.a_dag = a_dag;
InputStruct.kappa = kappa;
InputStruct.n_k = n_k;
InputStruct.hamilHandle = ...
    @ (t) calcHamiltonianTrans(n_hat, vPhi_hat, E_C, E_J);

%     @ (t) calcHamiltonianOsc(a, a_dag, w_r, hbar)
%     @ (t) calcHamiltonianDrive(a, a_dag, w_r, kappa, w_d, phi_d, ampHandle, t, hbar);
%     @ (t) calcHamiltonianTrans(n_hat, Vphi_hat, E_C, E_J);

%% Calculate Solution


rho_0 = reshape(rho, numel(rho), 1);
% funcHandle = @(t, y) masterEq(y, a, a_dag, w_r, kappa, n_k, hbar);
funcHandle = @(t, y) masterEq(t, y, InputStruct);
[t_out, y_out] = ode45(funcHandle, [0, 0.1], rho_0);



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