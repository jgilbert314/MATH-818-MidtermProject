function [ rho_dot ] = masterEq(t, rho_0, InputStruct)
% Summary:
%   This function defines an ODE for the state density matrix of a harmonic
%   oscillator connected to a transmission line.
% 
% Input:
%   t - time, dummy input for ode45
%   rho_0 - column vector representation of density matrix. Constructed 
%           from square matrix via:
%               rho_0 = reshape(rho, numel(rho), 1);
%   InputStruct - structure containing the following fields
%       a, a_dag - creation and ahnihilation operators (square matrices)
%                  a_dag = ctranspose(a);
%       kappa - damping coefficient
%       n_k - Constant determined by Bose-Einstein distribution
%                   From paper (p.6): n_k ~ 1.8
%       hamilHandle - a function handle defining a Hamiltonian: H(t)
% 
% Output:
%   rho_dot - time derivative of density matrix (column vector)

% Author: Jason Gilbert
% Date: March 15, 2022
% Version: V00
% Last Updated: March 16, 2022

% TODO:
%   ^ Rewrite to use InputStruct
%   ^ Allow variable Hamiltonian

%% Initialization

% Check that rho corresponds to square matrix
N = length(rho_0);
if ( mod(sqrt(N), 1) ~= 0 )
    error('Input vector must be perfect square in length');
end

a = InputStruct.a;
a_dag = InputStruct.a_dag;
kappa = InputStruct.kappa;
n_k = InputStruct.n_k;
hamilHandle = InputStruct.hamilHandle;


%% Calculation

% Construct square matrix
rho = reshape( rho_0, [sqrt(N), sqrt(N)] );


H = hamilHandle(t); % Oscillator Hamiltionian
Hcomm = H*rho - rho*H; % Commutator of H_s and rho

% Dissipative terms for master equation
diss1 = calcDissMat(a, rho);
diss2 = calcDissMat(a_dag, rho);

% Calculate Master equation
rho_dotMat = -1i*Hcomm + kappa*(n_k + 1)*diss1*rho + kappa*n_k*diss2*rho;

% Reshape to vector to work with ode solvers
rho_dot = reshape(rho_dotMat, [N, 1]);


end