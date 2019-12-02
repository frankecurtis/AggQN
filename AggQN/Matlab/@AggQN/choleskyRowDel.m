function R1 = choleskyRowDel(AQN,R,n)

% Author      : Baoyu ZHOU
% Description : Update Cholesky Factorizations by Row Deletion
% Input       : R       ~ initial Cholesky factor
%               n       ~ the index of column vector deleting from R'*R
% Output      : R1      ~ Cholesky factorization matrix after update

% size of initial Cholesky factor
s0 = size(R,1);

% divide initial Cholesky factor by parts
R11 = R(1:n-1,1:n-1);
R13 = R(1:n-1,n+1:s0);
R33 = R(n+1:s0,n+1:s0);
r23 = R(n,n+1:s0);

% Input matrix for Givens rotation
A = [r23;R33];

% Givens Rotation update
for i=1:s0-n
    alpha = A(i,i);
    beta = A(i+1,i);
    gamma = alpha / sqrt(alpha^2+beta^2);
    sigma = beta / sqrt(alpha^2+beta^2);
    T = [gamma sigma; -sigma gamma];
    A(i:i+1,i:s0-n) = T * A(i:i+1,i:s0-n);
    A(i+1,i) = 0;
end

% delete last row of zero
Rt33 = A(1:s0-n,:);

% update Cholesky factorization's result for row deletion operation
R1 = [R11 R13; zeros(s0-n,n-1) Rt33];

