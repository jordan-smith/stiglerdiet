function [x,y,obj] = solveLPRevised(A,b,c)
%
% This function executes phases I and II of the revised simplex method,
% returning an optimal solution if the LP is feasible with a bounded
% objective.
%
% On input:
% A: matrix of coefficients in the LP constraints
% b: requirement (column) vector
% c: cost (column) vector
%

[m n] = size(A);
Aext = [A eye(m)];
tempx = [zeros(n,1); b];
tempc = [zeros(1,n) ones(1,m)]'; 
tempB = find(tempx);
format short g;
disp('Starting Phase I...');

[newx newy newB flg tempt temps] = ComputePhase2(Aext,b,tempc,tempx,tempB);
disp('Phase I complete...');

if (flg == 1)
    disp('failure: an unknown problem occured');
end

%Case 2:
obj = tempc'*newx;
if (obj > 0)
    disp('LP has no feasible solution');
    disp('The optimal solution obtained in Phase I is');
    disp(newx);
    disp('objective value:');
    disp(obj);
    return;
end

%Case 1: obj == 0
disp('Starting Phase II...');
tempx = newx(1:n);
[x y flg t s] = ComputePhase2(A,b,c,tempx,newB);
disp('Phase II complete...');

if (flg == 1)    
    disp('LP is feasible with unbounded objective');
    disp('Displaying the set of feasible solutions x(L) that leads to -inf');
    B = find(x);
    for c = 1:n
        if (c == s)
            disp('lambda');
            continue;
        end
        index = find(B==c,1);
        if (x(c) > 0)
            if (t(index) < 0)
                nt = -t(index);
                disp([num2str(x(c)) ' + ' num2str(nt) '*lambda']);
                continue;
            end
            disp(num2str(x(c)));
            continue;
        end
        disp(num2str(0));
    end
    return;
end
obj = c'*x;
disp('LP is feasible with bounded objective');
disp('Displaying optimal primal solution:');
disp(x);
disp('Displaying optimal dual solution:');
disp(y);
disp(['optimal objective value = ' num2str(obj)]);




