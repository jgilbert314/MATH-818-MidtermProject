%% Initializaiton
clear;
close all;
DefineConstants; % Run initialization script

psi = [0:3]; % State vectors
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
n_hat = Q_hat/2/e; % Defined below Eq.(22)
vPhi_hat = 2*e/hbar*Phi_hat; % Phi_0 defined below Eq.(19)

% Ladder operators
[b, b_dag] = calcTransmonLadderOp(a, a_dag, Z_r, E_C, E_J, hbar, e);



% Density matrix
% psi_0 = zeros(1, N);
% psi_1 = psi_0;
% psi_1(1) = 1;
% psi_2 = psi_0;
% psi_2(2) = 1;
rho = kron(psi', psi); % Density matrix


%% Initialize Input Structure
model_flag = driven_oscillator;

w_d = w_r;
phi_d = pi*1.15;
ampHandle = @(t) 10/hbar;


g1 = 0.1*w_r;
g2 = g1;
E_C1 = E_C;
E_C2 = E_C1;
E_J1 = E_J;
E_J2 = E_J1;

[b1, b1_dag] = calcTransmonLadderOp(a, a_dag, Z_r, E_C1, E_J1, hbar, e);
[b2, b2_dag] = calcTransmonLadderOp(a, a_dag, Z_r, E_C2, E_J2, hbar, e);

InputStruct.a = a;
InputStruct.a_dag = a_dag;
InputStruct.kappa = kappa;
InputStruct.n_k = n_k;
if ( strcmp(model_flag, pure_oscillator) )
    disp('Running Pure Oscillator');
    InputStruct.hamilHandle = ...
        @ (t) calcHamiltonianOsc(a, a_dag, w_r, hbar);
elseif ( strcmp(model_flag, driven_oscillator) )
    disp('Running Driven Oscillator');
    InputStruct.hamilHandle = ...
        @ (t) calcHamiltonianDrive(a, a_dag, w_r, kappa, w_d, phi_d, ampHandle, t, hbar);
elseif ( strcmp(model_flag, single_transmon) )
    disp('Running Single Transmon');
    InputStruct.hamilHandle = ...
        @ (t) calcHamiltonianTrans(n_hat, vPhi_hat, E_C, E_J);
elseif ( strcmp(model_flag, rescoupled_transmon) )
    disp('Running Resonantly Coupled Transmon');
    [b1_eigvec, ~] = eigs(b1);
    psi_1 = b1_eigvec(:, 1);
    [b2_eigvec, ~] = eigs(b2);
    psi_2 = b2_eigvec(:, 1);
    rho = kron(psi_1', psi_2); % Density matrix
    
    InputStruct.hamilHandle = ...
        @ (t) calcOscCoupledTransmons(b1, b1_dag, b2, b2_dag, g1, g2, E_C1, E_J1, E_C2, E_J2, a, a_dag, w_r, hbar);
end


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


%%

[N_t, N_rho] = size(y_out);
N_d = round(sqrt(N_rho));
y_eig = zeros(N_d, N_t);

for itr = 1:N_t
   this_rho = reshape(y_out(itr, :), [N_d, N_d]);
   y_eig(:, itr) = sort( eig(this_rho) );
end

figure(4);
plot( t_out, abs(y_eig') );