function [X_minus,Qxx_minus] = ffun(X,Qxx,Phi,Qww)
% function to derive the predcited state vector and its VCM
%
%   INPUT:
%       X: state vector
%       Qxx: VCM of the states
%       Phi: transition matrix
%       Qww: VCM of the process noise
%
%   OUTPUT:
%       X_minus: predicted state vector
%       Qxx_minus: predicted VCM of the states
X_minus = Phi * X;                     % [4*4] * [4*1]
Qxx_minus = Phi * Qxx * Phi'+ Qww;     % [4*4] * [4*4] * [4*4] + [4*4]
    
    
end