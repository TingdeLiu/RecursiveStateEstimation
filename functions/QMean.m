function Qxx_minus = QMean(delta,W)
%   function to derive the predicted VCM of the states
%
%   INPUT:
%       delta: difference bewteen the Sigma points and the mean of the
%       Sigma points, which represents the state vector. dimention is 
%       [n x NoSig] with n: number of the unknown parameters (size of the 
%       state vector) and NoSig: number of the Sigma points in total
%       W: weight of the Sigma points with [NoSig x 1] dimension
%
%   OUTPUT:
%       Qxx_minus: predicted VCM of the states
Qxx_minus = zeros(size(delta,1),size(delta,1));
for i = 1:size(W,1)
    Qxx_i = W(i,1) * (delta(:,i) * delta(:,i)'); %[1*1]*[4*1]*[1*4]
    Qxx_minus = Qxx_minus + Qxx_i;
end

end