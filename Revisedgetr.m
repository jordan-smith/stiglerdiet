function r = Revisedgetr(n,s,B,T,t)
%
% find the index to kick out
%
% On input: 
% B: current basis index
% T: current Revised Tableau
% t: current pivot column
% s: s-th  column is to JOIN the index.
% n: number of unknowns.
%

% On output: 
%
% r: B(r)-th column is to LEAVE the index.
%    r < 0 indicates unbounded LP.
%

x   = zeros(n,1);
x(B)= T(1:end-1,1);
if (max(t)<n*eps)
    %Case 2a: no lower bound
    
    r = -1;
    return;
end
mask        = find(t>0);
posL        = x(B(mask))./t(mask);
[lamda, r]  = min(posL);
r           = mask(r);


