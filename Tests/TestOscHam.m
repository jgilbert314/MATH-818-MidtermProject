%% Initialization
clear;
close all;
DefineConstants;

save_flag = 0;
filename = fullfile(file_folder, 'OscTest');

%% Define constants
% Defines various constants used in the model
% Physical Constants

e = 1.602e-19;        % Charge of electron [C]
hbar = 1.0545718e-34; % Reduced Planck's constant [m^2*kg/s]

% Oscillator parameters
w_r = 2*pi*7e9;       % Cavity frequency
Z_r = 50;             % Transmission Line Impedance [Ohms]
n_k = 0;              % Constant from Bose-Einstein distribution
kappa = 0;            % Photon decay rate (damping factor)

%%
N = 5;
psi = [0:N]; % State vectors


% Calc ladder operators
[a, a_dag] = calcOscLadderOp(N);

% Oscillator Basis
[Phi_hat, Q_hat] = calcOscBasis(a, a_dag, Z_r, hbar);

% Density matrix
psi = psi/norm(psi); % Normalize so tr(rho^2) = 1
rho = kron(psi', psi); % Density matrix


%% Initialize Input Structure

InputStruct.a = a;
InputStruct.a_dag = a_dag;
InputStruct.kappa = kappa;
InputStruct.n_k = n_k;

disp('Running Pure Oscillator');
InputStruct.hamilHandle = ...
    @ (t) calcOscHamiltonian(a, a_dag, w_r, hbar);


%% Calculate Solution

rho_0 = reshape(rho, numel(rho), 1);
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