function h = hfun(TS1,TS2,X)
%   function to estimate the osbervations based on the most recent state
%   vector
%
%   INPUT:
%       TS1: position of the first total station (TS1)
%       TS2: position of the second total station (TS2)
%       X: state vector
%
%   OUTPUT:
%       h: estimated observation vector
d_TS1 = sqrt( (X(1,1)-TS1(1,1))^2 + (X(2,1)-TS1(2,1))^2 );
d_TS2 = sqrt( (X(1,1)-TS2(1,1))^2 + (X(2,1)-TS2(2,1))^2 );
alpha_TS1 = atan2(X(2,1)-TS1(2,1),X(1,1)-TS1(1,1));
alpha_TS2 = atan2(X(2,1)-TS2(2,1),X(1,1)-TS2(1,1));

h=[d_TS1,alpha_TS1,d_TS2,alpha_TS2];     
        
end