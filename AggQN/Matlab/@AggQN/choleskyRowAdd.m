function R1 = choleskyRowAdd(AQN,R,cv,n)

% Author      : Baoyu ZHOU
% Description : Update Cholesky Factorizations by Row Addition
% Input       : R       ~ initial Cholesky factor
%               cv      ~ column vector
%               n       ~ the index for v plugging into R'*R
% Output      : R1      ~ Cholesky factorization matrix after update


% size of initial Cholesky factor
s0 = size(R,1);

% divide initial Cholesky factor by parts
R11 = R(1:n-1,1:n-1);
R12 = R(1:n-1,n:s0);
R22 = R(n:s0,n:s0);

% divide row vector v by parts
a = cv(1:n-1);
b = cv(n);
c = cv(n+1:s0+1);

% calculate parts of new Cholesky factor
u = R11'\a;
v = sqrt(b-norm(u)*norm(u));
w = (c-R12'*u)/v;
w0 = w;

% set temporary variable
Rt22 = R22;

% get length of vector w
t = length(w);

% update Cholesky factorization
for k = 1:t
    a1 = sqrt(Rt22(k,k)^2 - w(k)^2);
    a2 = a1 / Rt22(k,k);
    a3 = w(k) / Rt22(k,k);
    Rt22(k,k) = a1;
    Rt22(k,k+1:t) = (Rt22(k,k+1:t) - a3*w(k+1:t,1)') / a2;
    w(k+1:t) = a2*w(k+1:t) - a3*Rt22(k,k+1:t)';
end

% update Cholesky factorization's result for row addition operation
R1 = [R11 u R12;zeros(1,n-1) v w0';zeros(s0-n+1,n-1) zeros(s0-n+1,1) Rt22];
