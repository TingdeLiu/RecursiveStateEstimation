function d = d_Vector(R,D,X,TS)
%   function to derive the "d" vector related to the constraints based on
%   the Dx = d constraint equation
%
%   INPUT:
%       R: constraint information
%       D: the derived D-Matrix
%       X: state vector
%       TS: position of the total station to which the contraint is related
%
%   OUTPUT:
%       d: d-vector [c x 1] dimension -> c: number of the constraints
         
syms x y 
constraint_equ = sqrt( (x-TS(1))^2 + (y-TS(2))^2 );

x=X(1,1); y=X(2,1);
d=R-eval(constraint_equ)+D*X;   
end