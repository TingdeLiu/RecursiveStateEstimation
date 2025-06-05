function X_minus = ffunStates(X,Phi)
% function to derive the predcited state vector
%
%   INPUT:
%       X: state vector
%       Phi: transition matrix
%
%   OUTPUT:
%       X_minus: predicted state vector
X_minus = Phi * X;    
end