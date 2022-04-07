%% Initializaiton
clear;
% close all;
DefineConstants; % Run initialization script

save_flag = 1;
filename_base = fullfile( file_folder, 'RefFig6' );

% Transmon Parameters
w_p = 2*pi*5e9;       % Plasma frequency sqrt(8*E_C*E_J)/hbar (p.11)
beta_vec = [2, 10, 50, 1e5]; % Ratio E_J/E_C
E_C = w_p*hbar/sqrt(8*beta); % Josephson energy
E_J = beta*E_C;              % Charging energy

w_q = calcQubitFreq(E_C, E_J, hbar); % Qubit frequency (p.10)
C_q = e^2/2/E_C;
L_q = (hbar/2/e)^2/E_J;
Z_r = sqrt(1/L_q/C_q);
% Z_r = hbar/e^2*sqrt(1/2/beta);
Z_r = 50; % Transmission Line Impedance [Ohms]

%% Basis Definition

N = 2;
psi = [0:N]; % State vectors


% Calc ladder operators
[a, a_dag] = calcOscLadderOp(N);
[b, b_dag] = calcTransLadderOp(a, a_dag, Z_r, E_C, E_J, hbar, e);


% Oscillator Basis
[Phi_hat, Q_hat] = calcOscBasis(a, a_dag, Z_r, hbar);
[n_hat, vPhi_hat] = calcTransBasis(Phi_hat, Q_hat, hbar, e);

%%
N_n = 1e2;
N_g = 1;
n_g = linspace(-N_g, N_g, N_n);

for itrB = 1:length(beta_vec)
    beta = beta_vec(itrB);
    E_J = beta*E_C;              % Charging energy

    eig_vals = zeros(N+1, N_n);
    cos_mat = ( expm(1i*vPhi_hat) + expm(-1i*vPhi_hat) )/2;
    for itr = 1:N_n
        H_hat = 4*E_C*(n_hat - n_g(itr))^2 - E_J*cos_mat;
        
        these_eig = sort( eig(H_hat) );
        eig_vals(:, itr) = these_eig - these_eig(1);
    end
    eig_vals = eig_vals/hbar/2/pi * 1e-9; % Convert to GHz
    
    % Add results to plot
    subplot(2, 2, itrB);
    plot(n_g, abs(eig_vals(:, :))', '.');
    xlabel('n_g');
    ylabel('Energy (GHz)');
    title(['E_J/E_C = ', num2str(beta)]);
    grid on;
    drawnow;
    
    
end

if (save_flag)
    print2pdf(gcf, [filename_base]);
end
