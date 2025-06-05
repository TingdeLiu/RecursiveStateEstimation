function D = D_Matrix(X,TS)
%   function to derive the "D" matrix related to the constraints based on
%   the Dx = d constraint equation
%
%   INPUT:
%       X: state vector
%       TS: position of the total station to which the constraint is
%       related
%
%   OUTPUT:
%       D: D-Matrix with [c x u] dimension -> c: number of the constraints
%       and u: number of the unknowns

syms x y 
constraint_equ = sqrt( (x-TS(1))^2 + (y-TS(2))^2 );
dx_constraint_equ = diff(constraint_equ,x);
dy_constraint_equ = diff(constraint_equ,y);

x=X(1,1); y=X(2,1);

D=[eval(dx_constraint_equ),eval(dy_constraint_equ),0,0];
end