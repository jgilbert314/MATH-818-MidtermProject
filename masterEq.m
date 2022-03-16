function [ rho_dot ] = masterEq(rho_0, a, a_dag, w_r, kappa, n_k, hbar)
% Summary:
%   This function defines an ODE for the state density matrix of a harmonic
%   oscillator connected to a transmission line.
% 
% Input:
%   rho_0 - column vector representation of density matrix. Constructed 
%           from square matrix via:
%               rho_0 = reshape(rho, numel(rho), 1);
%   a, a_dag - creation and ahnihilation operators (square matrices)
%               a_dag = ctranspose(a);
%   w_r - resonant frequency of microwave cavity
%               From paper (p.5): w_r ~ 2*pi * 8 GHz
%   kappa - damping coefficient
%   n_k - Constant determined by Bose-Einstein distribution
%               From paper (p.6): n_k ~ 1.8
%   hbar - reduced Planck's constant
%           hbar ~ 1.0545718 × 10-34 [m^2*kg/s]
% 
% Output:
%   rho_dot - time derivative of density matrix (column vector)

% Author: Jason Gilbert
% Date: March 15, 2022
% Version: V00
% Last Updated: N/A


% Check that rho corresponds to square matrix
N = length(rho_0);
if ( mod(sqrt(N), 1) ~= 0 )
    error('Input vector must be perfect square in length');
end

% Construct square matrix
rho = reshape( rho_0, [sqrt(N), sqrt(N)] );


H_s = calcHamiltonianOsc(a, a_dag, w_r, hbar); % Oscillator Hamiltionian
Hcomm = H_s*rho - rho*H_s; % Commutator of H_s and rho

% Dissipative terms for master equation
diss1 = calcDissMat(a, rho);
diss2 = calcDissMat(a_dag, rho);

% Calculate Master equation
rho_dotMat = -1i*Hcomm + kappa*(n_k + 1)*diss1*rho + kappa*n_k*diss2*rho;

% Reshape to vector to work with ode solvers
rho_dot = reshape(rho_dotMat, [N, 1]);


end