% Copyright (C) 2021 Albert Berahas, Frank E. Curtis, Baoyu Zhou
%
% All Rights Reserved.
%
% Authors: Albert Berahas, Frank E. Curtis, Baoyu Zhou
%
% Method definition for AggQN class

% Compute 'b', "upper" part of 'A', 'Omega', 'omega', 'rhs' values, and auxiliary values
function computeAggregationValuesInitial(AQN)

% Construct HS_j
AQN.HS_j = AQN.computeInnerHessianProduct(AQN.S(:,AQN.j:end));

% Construct SHS_j
AQN.SHS_j = AQN.S(:,AQN.j:end)'*AQN.HS_j;

% Compute factorization of Q'*H*Q
[L_SHS,~] = AQN.choleskyPerturb(AQN.S(:,AQN.j:end)'*AQN.HS_j);

% Compute Q^{-1}
AQN.invM_j = L_SHS\eye(size(AQN.Y,2)-AQN.j+1);
AQN.invM_j = L_SHS'\AQN.invM_j;

% Compute b_j = -rho_j * (S(:,j:end)'*Y(:,j:end-1) - R(:,j:end-1))' * tau_j
AQN.b_j = -AQN.rho_j * tril(AQN.SY(AQN.j:end,AQN.j:end-1),-1)' * AQN.tau_j;

% Initialize matrix
AQN.MA_j = zeros(size(AQN.Y,2)-AQN.j+1,size(AQN.Y,2)-AQN.j);

% Compute "top" part of aggregation matrix M*A (eq. (3.15))
for i = 1:size(AQN.Y,2)-AQN.j
  AQN.MA_j(1:i,i) = -AQN.b_j(i)*AQN.S(:,AQN.j:AQN.j+i-1)'*AQN.y_j;
end

% Call unit test for "top" part computation (eq. (3.13))
AQN.runUnitTest(1);

% Compute Omega
AQN.Omega_j = AQN.S(:,AQN.j:end)'*AQN.y_j*AQN.b_j' + AQN.S(:,AQN.j:end)'*AQN.Y(:,AQN.j:end-1) - triu(AQN.SY(AQN.j:end,AQN.j:end-1),0);

% Compute omega
AQN.omega_j = AQN.b_j/sqrt(AQN.rho_j);

% Compute right-hand side
AQN.rhs_j = AQN.omega_j * AQN.omega_j' + AQN.Omega_j' * AQN.invM_j * AQN.Omega_j;

end
