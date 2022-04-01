%% Define constants
% Defines various constants used in the model

% Physical Constants
e = 1.602e-19;        % Charge of electron [C]
hbar = 1.0545718e-34; % Reduced Planck's constant [m^2*kg/s]

% Oscillator parameters
w_r = 2*pi*5e9;       % Cavity frequency
Z_r = 50;             % Transmission Line Impedance [Ohms]
n_k = 1.8;            % Constant from Bose-Einstein distribution
kappa = 10;            % Photon decay rate

% Transmon Parameters
w_p = 2*pi*5e9;       % Plasma frequency sqrt(8*E_C*E_J)/hbar (p.11)
beta = 1e6;             % Ratio E_J/E_C
E_C = w_p*hbar/sqrt(8*beta); % Josephson energy
E_J = beta*E_C;              % Charging energy


%% Model flags
% Strings specifying different model choices

pure_oscillator = 'Pure Oscillator';
driven_oscillator = 'Driven Oscillator';
single_transmon = 'Single Transmon';
rescoupled_transmon = 'Resonator Coupled Transmon';

