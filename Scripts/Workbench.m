%% Initializaiton
clear;
close all;
DefineConstants; % Run initialization script

%% Model flags
% Strings specifying different model choices

pure_oscillator = 'Pure Oscillator';
driven_oscillator = 'Driven Oscillator';
single_transmon = 'Single Transmon';
rescoupled_transmon = 'Resonator Coupled Transmon';

%%
N = 5;
psi = [0:N]; % State vectors


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
% psi_0 = zeros(1, N);
% psi_1 = psi_0;
% psi_1(1) = 1;
% psi_2 = psi_0;
% psi_2(2) = 1;
rho = kron(psi', psi); % Density matrix


%% Initialize Input Structure
model_flag = rescoupled_transmon;

w_d = w_r;
phi_d = pi*1.15;
ampHandle = @(t) 10/hbar;


g1 = 0.1*w_r;
g2 = g1;
E_C1 = E_C;
E_C2 = E_C1;
E_J1 = E_J;
E_J2 = E_J1;

[b1, b1_dag] = calcTransLadderOp(a, a_dag, Z_r, E_C1, E_J1, hbar, e);
[b2, b2_dag] = calcTransLadderOp(a, a_dag, Z_r, E_C2, E_J2, hbar, e);

InputStruct.a = a;
InputStruct.a_dag = a_dag;
InputStruct.kappa = kappa;
InputStruct.n_k = n_k;
if ( strcmp(model_flag, pure_oscillator) )
    disp('Running Pure Oscillator');
    InputStruct.hamilHandle = ...
        @ (t) calcOscHamiltonian(a, a_dag, w_r, hbar);
elseif ( strcmp(model_flag, driven_oscillator) )
    disp('Running Driven Oscillator');
    InputStruct.hamilHandle = ...
        @ (t) calcDriveHamiltonian(a, a_dag, w_r, kappa, w_d, phi_d, ampHandle, t, hbar);
elseif ( strcmp(model_flag, single_transmon) )
    disp('Running Single Transmon');
    InputStruct.hamilHandle = ...
        @ (t) calcTransHamiltonian(n_hat, vPhi_hat, E_C, E_J);
elseif ( strcmp(model_flag, rescoupled_transmon) )
    disp('Running Resonantly Coupled Transmon');
    % Ladder operators for master equation
    InputStruct.a = b;
    InputStruct.a_dag = b_dag;
    % Eigenstates of first transmon
    [b1_eigvec, ~] = eigs(b1); 
    psi_1 = b1_eigvec(:, 1);
    % Eigenstates of second transmon
    [b2_eigvec, ~] = eigs(b2);
    psi_2 = b2_eigvec(:, 1);
    % Density matrix
    rho = kron(psi_1', psi_2);
    
    InputStruct.hamilHandle = ...
        @ (t) calcOscCoupledTransHamiltonian(g1, g2, E_C1, E_J1, E_C2, E_J2, a, a_dag, w_r, Z_r, hbar, e);
end


%% Calculate Solution

rho_0 = reshape(rho, numel(rho), 1);
% funcHandle = @(t, y) masterEq(y, a, a_dag, w_r, kappa, n_k, hbar);
funcHandle = @(t, y) masterEq(t, y, InputStruct);
[t_out, y_out] = ode45(funcHandle, [0, 10], rho_0);



%% Plot Solution

diag_inds = 1 + [0:N-1]*(N+1);
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



%%

test = eig(vPhi_hat);
test2 = this_rho*test;
test2 = test*test2';

imagesc( abs(test2) );