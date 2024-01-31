function [W,Z] = LRSA(X,lambda,beta) 

% min_{Z,B,J,W}
% ||W||_k+lambda*||X-XZ||_F^2+beta||J||_1+0.5*mu*(||Z-B||_F^2+||Z-J||_F^2)+tr(Y1'*(Z-J)+Y2'*(Z-B))
% s.t. W=(|B|+|B'|)/2 
% In order to make the problem solvable, variable A is added,
% and the problem is transformed to   
%min_{Z,B,J,W,A}
% tr((Diag(W1)-W)*A)+lambda*||X-XZ||_F^2+beta*||J||_1+0.5*mu*(||Z-B||_F^2+||Z-J||_F^2)+tr(Y1'*(Z-J)+Y2'*(Z-B))
% s.t. W=(|B|+|B'|)/2 
% 0<=A<=I, Tr(A)=n.

if nargin < 4
    display = 0;
end
tol = 1e-4;
n = size(X,2);
rho = 1.1;
max_mu = 1e6;
mu = 1e-6;
maxIter = 1e6; 
one = ones(n,1);
XtX = X'*X;
I = eye(n);

Z = zeros(n);
J = Z;
B = Z;
A = Z;

Y1 = zeros(n,n);
Y2 = zeros(n,n);

iter = 0;
DEBUG = 1;
if DEBUG
disp(['initial,rank=' num2str(rank(Z))]);
end
while iter < maxIter
    iter = iter + 1;       
    % update Z
    invXtXI = I/(2*lambda*XtX+2*mu*I);
    Z = invXtXI*(2*lambda*XtX-Y1-Y2+mu*J+mu*B);
    
    % update J
    M = Z+Y1/mu;
    temp1=beta/mu;
    J = max(M-temp1,0)+min(M+temp1,0);
    
    % update B
    N = Z+Y2/mu;
    P = diag(A)*one'-A;
    V = (P+P')/mu;
    B = solve_b(N,V,1);
    
    W = (abs(B)+abs(B'))/2;
    L = diag(W*one)-W;    
    
    % update A
    [V1, D] = eig(L);
    D = diag(D);
    [~, ind] = sort(D);    
    A = V1(:,ind(1:n))*V1(:,ind(1:n))';
    
    leq1 = Z-J;
    leq2 = Z-B;
    stopC = max(max(max(abs(leq1))),max(max(abs(leq2)))); 
    
    if DEBUG
    if iter==1 || mod(iter,50)==0 || stopC<tol
        disp(['iter ' num2str(iter) ',mu=' num2str(mu,'%2.1e') ...
            ',rank=' num2str(rank(Z,1e-3*norm(Z,2))) ',stopALM=' num2str(stopC,'%2.3e')]);
    end
    end
    
    if stopC < tol 
        break;
    else
        Y1 = Y1 + mu*leq1;
        Y2 = Y2 + mu*leq2;
        mu = min(max_mu,mu*rho);
    end
end
end


function [B]=solve_b(N,V,gamma)
[a,b]=size(N);
B=zeros(a,b);
for i=1:a
    for j=1:b
        if V(i,j)>=0
        B(i,j)=max(N(i,j)-gamma*V(i,j)/2,0)+min(N(i,j)+gamma*V(i,j)/2,0);
        else
            if N(i,j)>=0
                B(i,j)=N(i,j)-gamma*V(i,j)/2;
            else
                B(i,j)=N(i,j)+gamma*V(i,j)/2;
            end
        end
    end  
end
end
