function [SigmaP,W] = UT_transform(NoSig,X,Qxx,gamma,n)
%   function to derive the deterministic Sigma points needed in the
%   Unscented Kalman Filter (UKF) and their corresponding weights
%
%   INPUT:
%       NoSig: number of the Sigma points that should be generated
%       X: state vector
%       Qxx: VCM of the states
%       gamma: scale factor in UKF
%       n: number of the unknown parameters (size of the state vector)
%
%   OUTPUT:
%       SigmaP: the generated Sigma points with [n x NoSig] dimension
%       W: weight of the generated Sigma points with [NoSig x 1] dimension

SigmaP = zeros(n,NoSig);
W = zeros(NoSig,1);
SigmaP(:,1) = X;      % x0=X
lambda = gamma ^ 2  - n;
W(1,1) = lambda / ( n + lambda );           % w0=0
deviation = (gamma * chol(Qxx))';  % where (.) i denotes the ith column of the matrix square root

% 1..nx
for i = 1:n   
SigmaP(:,i+1) = X + deviation(:,i);
W(i+1,1) = 1 / ( 2 * (n + lambda) );
end

% nx..2nx
for i = 1:n  
SigmaP(:,i+1+n) = X - deviation(:,i);
W(i+1+n,1) = 1 / ( 2 * (n + lambda) );
end


end
