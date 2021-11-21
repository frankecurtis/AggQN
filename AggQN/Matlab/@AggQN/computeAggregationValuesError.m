% Copyright (C) 2021 Albert Berahas, Frank E. Curtis, Baoyu Zhou
%
% All Rights Reserved.
%
% Authors: Albert Berahas, Frank E. Curtis, Baoyu Zhou
%
% Method definition for AggQN class

% Evaluates aggregation error
function [errorMaxAbs,errorFroSqr,v1,v2] = computeAggregationValuesError(AQN,A_j)

% Compute Ytilde
Ytilde = AQN.Y;
Ytilde(:,AQN.j:end-1) = AQN.Y(:,AQN.j:end-1) + AQN.computeInnerHessianProduct(AQN.S(:,AQN.j:end)*A_j) + AQN.y_j*AQN.b_j';

% Compute "full" b
b = [zeros(AQN.j-1,1); AQN.b_j; 0];

% Compute "full" A
A = [A_j zeros(size(AQN.S(:,AQN.j:end),2),1)];

% Construct matrices
matrix1 = triu(AQN.SY(AQN.j:end,AQN.j:end)) - triu(AQN.S(:,AQN.j:end)'*Ytilde(:,AQN.j:end));
matrix2 = b(AQN.j:end) + AQN.rho_j*tril(AQN.SY(AQN.j:end,AQN.j:end),-1)'*AQN.tau_j;
matrix3 = (Ytilde(:,AQN.j:end) - AQN.Y(:,AQN.j:end))'*AQN.computeInnerInverseHessianProduct(Ytilde(:,AQN.j:end)-AQN.Y(:,AQN.j:end)) - ...
  (1+AQN.rho_j*AQN.y_j'*AQN.computeInnerInverseHessianProduct(AQN.y_j))/AQN.rho_j * (b(AQN.j:end)*b(AQN.j:end)') + ...
  A'*tril(AQN.SY(AQN.j:end,AQN.j:end),-1) + tril(AQN.SY(AQN.j:end,AQN.j:end),-1)'*A;

% Check if trying Newton
if AQN.tryNewton

  % Set vectors for Newton's method
  v1 = [];
  for i = 1:size(AQN.Y,2)-AQN.j
    v1 = [v1 ; AQN.SHS_j(1:i,:)*A_j(:,i) + AQN.b_j(i)*AQN.S(:,1:i)'*AQN.y_j];
  end

  % Second: vectorize the lower part of the difference matrix
  M  = A_j'*AQN.SHS_j*A_j + AQN.Omega_j'*A_j + A_j'*AQN.Omega_j - AQN.omega_j*AQN.omega_j';
  v2 = zeros(size(M,1)*(size(M,1) + 1)/2,1);
  j  = 1;
  for i = 1:size(M,1)
    v2(j:j+i-1) = M(i,1:i)'; j = j+i;
  end

else

  % Empty return values
  v1 = []; v2 = [];

end

% Set errors
error1 = max(max(abs(matrix1)));
error2 =     max(abs(matrix2)) ;
error3 = max(max(abs(matrix3)));

% Set overall errors
errorMaxAbs = max([error1,error2,error3]);
errorFroSqr = (1/2)*(norm(matrix1,'fro')^2 + norm(matrix3,'fro')^2);

% Print errors
if AQN.verbosity >= 1
  fprintf('AggQN: Error in satisfaction of R = Rtilde             : %e\n',error1);
  fprintf('AggQN: Error in satisfaction of b definition           : %e\n',error2);
  fprintf('AggQN: Error in satisfaction of quadratic              : %e\n',error3);
  fprintf('AggQN: Error in satisfaction of all equations (MaxAbs) : %e\n',errorMaxAbs);
  fprintf('AggQN: Error in satisfaction of all equations (FroSqr) : %e\n',errorFroSqr);
end

end
