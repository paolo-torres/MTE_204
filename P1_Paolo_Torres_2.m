%{
Circuit 4

Mesh Equations:
6*i1 - 4*i2 - 2*i3  = 20
4*i1 - 18*i2 + 8*i3 =  0
2*i1 + 8*i2 - 15*i3 =  0

System of Equations (re-arranged, not diagonally dominant):
| 6  -4 -2  |  | i1 |  | 20  |
| 4 -18  8  |  | i2 |  |  0  |
| 2  -8 -15 |  | i3 |  |  0  |
%}

% clear output
clear;
clc;

% Gauss Seidel
i = 1;

% initial guesses
i2(i) = 0;
i3(i) = 0;

% initialize all errors to 100 to initiate while loop
error_i1(i) = 100;
error_i2(i) = 100;
error_i3(i) = 100;

% loop until the approximate relative error is less than or equal to 5%
while(error_i1(i) >= 0.05 || error_i2(i) >= 0.05 || error_i3(i) >= 0.05)
    % calculate next values for each variable
    i1(i+1) = (20 + 4.*i2(i) + 2.*i3(i))./6;
    i2(i+1) = (-4.*i1(i+1) - 8.*i3(i))./(-18);
    i3(i+1) = (-2.*i1(i+1) - 8.*i2(i+1))./(-15);

    % calculate next errors for each variable
    error_i1(i+1) = abs((i1(i+1)-i1(i))./i1(i+1)).*100;
    error_i2(i+1) = abs((i2(i+1)-i2(i))./i2(i+1)).*100;
    error_i3(i+1) = abs((i3(i+1)-i3(i))./i3(i+1)).*100;
    
    % increment number of iterations
    i = i+1;
end

% Equation Solver
syms x1 x2 x3
% each equation of the matrix
equations = [6.*x1-4.*x2-2.*x3==20, 4.*x1-18.*x2+8.*x3==0, 2.*x1+8.*x2-15.*x3==0];
% all variables in the matrix
variables = [x1 x2 x3];
% built-in equation solver
S = solve(equations, variables);

% Matrix Inversion
A = [6 -4 -2; 4 -18 8; 2 8 -15];
% inverse of matrix
A_inv = inv(A);
B = [20; 0; 0];
% [A][x]=B --> [x]=[A]^-1[B]
T = A_inv*B;

% display all variable values, errors, solver, and inversion
disp('Gauss Seidel');
disp(i1(i));
disp(i2(i));
disp(i3(i));
disp('Errors');
disp(error_i1(i));
disp(error_i2(i));
disp(error_i3(i));
disp('Iterations');
disp(i - 1);
disp('Equation Solver');
disp(double(S.x1));
disp(double(S.x2));
disp(double(S.x3));
disp('Matrix Inversion');
disp(T);
