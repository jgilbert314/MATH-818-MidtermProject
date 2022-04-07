function [ w_q ] = calcQubitFreq(E_C, E_J, hbar)
% Definition for qubit frequency found below Eq.27 (p.10)
% 
% Input:
%   E_C, E_J - charging, Josephson energies
%   hbar - reduced Planck's constant
% 
% Output:
%   w_q - qubit frequency [rad/s]

w_q = ( sqrt(8*E_C*E_J) - E_C )/hbar;

end