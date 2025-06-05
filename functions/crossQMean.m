function Q_LP = crossQMean(deltaX,deltaL,W)
%   function to derrive the cross VCM 
%
%   INPUT:
%       deltaX: difference bewteen the Sigma points and the mean of the
%       Sigma points, which represents the state vector. dimention is 
%       [n x NoSig] with n: number of the unknown parameters (size of the 
%       state vector) and NoSig: number of the Sigma points in total
%       deltaL: difference between the estimated observations based on the
%       Sigma points and the mean of the estimated observations by the 
%       Sigma points, which represents the observation vector
%       W: weight of the Sigma points with [NoSig x 1] dimension
Q_LP = zeros(size(deltaX,1),size(deltaX,1));
for i = 1:size(W,1)
    Q_LP_i = W(i,1) * (deltaX(:,i) * deltaL(:,i)');
    Q_LP = Q_LP + Q_LP_i;
end

end