%%
clear;
close all;
DefineConstants; % Run initialization script

% Basis size
N = 10;
Z_r = 20;

figure(1);
RefFig30;




% Basis size
N = 6;
Z_r = 8;

figure(2);
RefFig30;





%% Check Eigenvalues

Z_r_vec = linspace(1, 100, 1e2);

for itr = 1:length(Z_r_vec)
    Z_r = Z_r_vec(itr);
    RefFig30;

    % Plot Results
    tran_freq = w_q2/2/pi * 1e-9; % Convert to GHz
    eig_freq = eig_arr/hbar/(2*pi)^2 * 1e-9; % Convert to GHz
    plot(tran_freq, eig_freq(1:5, :)');
    xlabel('Second Qubit Frequency [GHz]');
    ylabel('Eigenfrequency');
    % xlim([10, 20]);
    ylim([0, 20]);
    grid on;
    title(num2str(Z_r));
    drawnow;
end