function A = design_matrix(X,TS1,TS2)
%   function to derive the design matrix
%   
%   INPUT:
%       X: state vector
%       TS1: position of the first total station (TS1)
%       TS2: position of the second total station (TS2)
%
%   OUTPUT:
%       A: design matrix with [n x u] dimension -> n: number of the
%       observations and u: number of the unknowns

% h-function
syms x y d_TS1 d_TS2 alpha_TS1 alpha_TS2
d_TS1 = sqrt( (x-TS1(1,1))^2 + (y-TS1(2,1))^2 );
d_TS2 = sqrt( (x-TS2(1,1))^2 + (y-TS2(2,1))^2 );
alpha_TS1 = atan2(y-TS1(2,1),x-TS1(1,1));
alpha_TS2 = atan2(y-TS2(2,1),x-TS2(1,1));

d_d_TS1_x=diff(d_TS1,x);
d_d_TS1_y=diff(d_TS1,y);
d_alpha_TS1_x=diff(alpha_TS1,x);
d_alpha_TS1_y=diff(alpha_TS1,y);

% partial derivative of h-function
d_d_TS2_x=diff(d_TS2,x);
d_d_TS2_y=diff(d_TS2,y);
d_alpha_TS2_x=diff(alpha_TS2,x);
d_alpha_TS2_y=diff(alpha_TS2,y);

% Assign value to the function
x=X(1,1); y=X(2,1);

A=[eval(d_d_TS1_x),eval(d_d_TS1_y),0,0;...
    eval(d_alpha_TS1_x),eval(d_alpha_TS1_y),0,0;...
    eval(d_d_TS2_x),eval(d_d_TS2_y),0,0;...
    eval(d_alpha_TS2_x),eval(d_alpha_TS2_y),0,0];
    
         
end