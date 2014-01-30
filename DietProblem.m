% 
% Computes the optimal solution to George Stigler's Diet Problem
% as developed in his paper.
%
% written by Jordan Smith
%


% the matrix in StiglerMatrix.xlsx' is the transpose of the diet LP matrix
At = xlsread('StiglerMatrix.xlsx');
% Introduce slack variables to allow surpluses of nutrients
A = [At' -eye(9)];
b = [3000 70 .8 12 5000 1.8 2.7 18 75]';
c = [ones(1,77) zeros(1,9)]';

[x y obj] = solveLPRevised(A,b,c);
