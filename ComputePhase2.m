function [x,y,B,flg,t,s] = ComputePhase2(A,b,c,initx,initB)
%
% This function executes Phase II of the revised simplex method.
%
% On input:
% A: matrix of coefficients in the LP constraints
% b: requirement (column) vector
% c: cost (column) vector
% initx: initial basic feasible solution
% initB: intial Basis for the first phase II iteration
%

[m n] = size(A);
% Construct the initial basis and tableau
x = initx;
B = initB;
T = A(:,B)\[b eye(m)];
y = T(:,2:end)'*c(B);
T = [T;[c'*x,y']];
flg = 0;
%
% Starting Simplex Method
%
f = ['Starting Phase II Simplex Iteration... '];
format short g;
disp(f);
disp('Initial Basis is');
disp(B');
disp('Displaying Initial solution x, c-A^T*y and their componentwise product');
disp([x c-A'*y x.*(c-A'*y)]);
simplex = 1;
ITER = 0;
obj = c'*x;
pause(2);

while (simplex == 1)
%
% determine the next s and r values.
%
   y = T(end,2:end)';
   [zmin,s] = min(c-A'*y); 
%
% check for convergence.
%
   if (abs(zmin) < 1e-14)
       %Case 1: y is feasible for the dual
       disp('Simplex Method has converged');
       simplex = 0;
       x = zeros(n,1);
       obj = c'*x;
       x(B) = T(1:end-1,1);
       disp('Displaying Optimal Basis');
       disp(B');
       disp('Displaying Optimal solution x, c-A^T*y and their componentwise product');
       disp([x c-A'*y x.*(c-A'*y)]);
       continue;
   end

   t = T(1:end-1,2:end)*A(:,s);
   r = Revisedgetr(n,s,B,T,t);
   if (r < 1)
       %Case 2a
       disp('LP has no lower bound');
       simplex = 0;
       flg = 1;
       continue;
   end
   x   = zeros(n,1);
   x(B)= T(1:end-1,1);
   ITER = ITER + 1;
   f = ['Iteration ', num2str(ITER), ' Obj ', num2str(c'*x), '. Smallest component in c-A^T*y: ', ... 
         num2str(zmin), ' at s =', num2str(s), '. Component r = ', num2str(r), ' out of basis'];
   disp(f);
   obj1 = c'*x;
%
% update the revised simplex tableau.
%
   [T,B1,flg]=RevisedSimplexTableau(B,r,s,t,zmin,T);      
   if (flg == 1)
       disp('LP is degenerate');
       simplex = 0;
       continue;
   end
   B   = B1;
   obj = obj1;
   disp('Current Basis is');
   disp(B');
   pause(1);
end




