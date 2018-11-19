%{
Circuit 3

Mesh Equations:
20*i1 - 15*i2         =  80
15*i1 - 50*i2 + 25*i3 =   0
        25*i2 - 45*i3 =  50

System of Equations (diagonally dominant):
| 20 -15   0 |  | i1 |  | 80 |
| 15 -50  25 |  | i2 |  |  0 |
|  0  25 -45 |  | i3 |  | 50 |
%}

% clear output
clear;
clc;

% Gauss Seidel
i = 1;

% initial guesses
i1(i) = 0;
i2(i) = 0;
i3(i) = 0;

% initialize all errors to 100 to initiate while loop
error_i1(i) = 100;
error_i2(i) = 100;
error_i3(i) = 100;

% loop until the approximate relative error is less than or equal to 5%
while(error_i1(i) >= 0.05 || error_i2(i) >= 0.05 || error_i3(i) >= 0.05)
    % calculate next values for each variable
    i1(i+1) = (80 + 15.*i2(i) + 2.*i3(i))./20;
    i2(i+1) = (-15.*i1(i+1) - 25.*i3(i))./(-50);
    i3(i+1) = (50 + -25.*i2(i+1))./(-45);

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
equations = [20.*x1-15.*x2==80, 15.*x1-50.*x2+25.*x3==0, 25.*x2-45.*x3==50];
% all variables in the matrix
variables = [x1 x2 x3];
% built-in equation solver
S = solve(equations, variables);

% Matrix Inversion
A = [20 -15 0; 15 -50 25; 0 25 -45];
% inverse of matrix
A_inv = inv(A);
B = [80; 0; 50];
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
