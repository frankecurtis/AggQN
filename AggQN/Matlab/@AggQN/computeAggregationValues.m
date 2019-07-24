% Copyright (C) 2019 Albert Berahas, Frank E. Curtis, Baoyu Zhou
%
% All Rights Reserved.
%
% Authors: Albert Berahas, Frank E. Curtis, Baoyu Zhou
%
% Method definition for AggQN class

% Compute aggregation values 'A' and 'b'
function computeAggregationValues(AQN)

% Construct H_{1:j-1} as AQN_H.hessian
H_matr = AQN.initHv(eye(AQN.n));
AQN_H  = AggQN('H',H_matr);
AQN_H.setVerbosity(0);
for i = 1:AQN.j-1
  AQN_H.addPair(AQN.S(:,i),AQN.Y(:,i));
end

% Construct HS_j
AQN.HS_j = AQN_H.hessian*AQN.S(:,AQN.j:end);

% Compute factorization of S'*H*S
[AQN.L_SHS,~] = AQN.choleskyPerturb(AQN.S(:,AQN.j:end)'*AQN.HS_j);

% Compute Pre-Condition term C
[temp1,Lambda,~] = svd(AQN.S(:,AQN.j:end)'*AQN.S(:,AQN.j:end));
AQN.C = temp1*(sqrt(Lambda)\eye(size(AQN.S,2)-AQN.j+1));

% Compute Q
AQN.Q = AQN.S(:,AQN.j:end) * AQN.C;

% Construct HSC_j
AQN.HQ_j = AQN_H.hessian * AQN.Q;

% Compute factorization of Q'*H*Q
[L_QHQ,~] = AQN.choleskyPerturb(AQN.Q'*AQN.HQ_j);


% Compute Q^{-1}
invM = L_QHQ\eye(size(AQN.Y,2)-AQN.j+1);
invM = L_QHQ'\invM;

%%%%%%%%%%%%%
% Compute b %
%%%%%%%%%%%%%

% Compute b_j = -rho_j * (S(:,j:end)'*Y(:,j:end-1) - R(:,j:end-1))' * tau_j
AQN.b_j = -AQN.rho_j * tril(AQN.SY(AQN.j:end,AQN.j:end-1),-1)' * AQN.tau_j;

%%%%%%%%%%%%%%%%%%%%%%%%
% Compute "top" of M*A %
%%%%%%%%%%%%%%%%%%%%%%%%

% Initialize matrix
AQN.MA_j = zeros(size(AQN.Y,2)-AQN.j+1,size(AQN.Y,2)-AQN.j);

% Compute "top" part of aggregation matrix M*A (eq. (3.15))
for i = 1:size(AQN.Y,2)-AQN.j
  AQN.MA_j(1:i,i) = -AQN.b_j(i)*AQN.S(:,AQN.j:AQN.j+i-1)'*AQN.y_j;
end

% Call unit test for "top" part computation (eq. (3.13))
AQN.runUnitTest(1,invM);

%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Compute "bottom" of Q*A %
%%%%%%%%%%%%%%%%%%%%%%%%%%%

% Compute Omega
Omega = AQN.Q'*AQN.y_j*AQN.b_j' + AQN.Q'*AQN.Y(:,AQN.j:end-1) - AQN.C'* triu(AQN.SY(AQN.j:end,AQN.j:end-1),0);

% Compute omega
omega = AQN.b_j/sqrt(AQN.rho_j);

% Compute right-hand side
rhs = omega * omega' + Omega' * invM * Omega;

% Solve for the first quadratic equation (eq. (3.21c)), solved directly
% without computing a_{m-1,2}^* + lambda_{m-1}
invMC = AQN.C*invM*AQN.C';
AQN.MA_j(size(AQN.Y,2)-AQN.j+1,size(AQN.Y,2)-AQN.j) = sqrt(rhs(size(AQN.Y,2)-AQN.j,size(AQN.Y,2)-AQN.j)/invMC(size(AQN.Y,2)-AQN.j+1,size(AQN.Y,2)-AQN.j+1)) - AQN.S(:,end)'*(AQN.b_j(size(AQN.Y,2)-AQN.j)*AQN.y_j + AQN.Y(:,size(AQN.Y,2)-1));

% Call unit test for first quadratic
AQN.runUnitTest(2,invMC,rhs);

% Initialize phi
phi = zeros(size(AQN.Y,2)-AQN.j+1,size(AQN.Y,2)-AQN.j-1);

% Initialize values for "cutting off" beta
good_ind = [];

% Loop through remaining levels in reverse order
for level = size(AQN.Y,2)-AQN.j-1:-1:1
  
  % Compute phi
  phi(level+2:size(AQN.Y,2)-AQN.j+1,level) = AQN.MA_j(level+2:size(AQN.Y,2)-AQN.j+1,level+1) + AQN.S(:,level+AQN.j+1:end)'*(AQN.b_j(level+1)*AQN.y_j + AQN.Y(:,level+AQN.j));
  
  % Compute beta
  good_ind = [level; good_ind];
  if cond(phi(:,good_ind)'*phi(:,good_ind)) <= AQN.cond_tol_beta
    beta_sub = (phi(:,good_ind)'*invMC*phi(:,good_ind))\rhs(good_ind+1,level);
    beta     = zeros(size(AQN.Y,2)-AQN.j-level,1);
    for i = 1:length(good_ind)
      beta(good_ind(i)+1-level) = beta_sub(i);
    end
  else
    good_ind = good_ind(2:end);
    beta_sub = (phi(:,good_ind)'*invMC*phi(:,good_ind))\rhs(good_ind+1,level);
    beta     = zeros(size(AQN.Y,2)-AQN.j-level,1);
    for i = 1:length(good_ind)
      beta(good_ind(i)+1-level) = beta_sub(i);
    end
  end

  % Compute combination of phi
  phi_comb = phi(:,level:size(AQN.Y,2)-AQN.j-1)*beta;

  % Call unit test for phi combination
  AQN.runUnitTest(3,invMC,rhs,level,phi,phi_comb);
  
  % Compute null space for affine equations
  N = null(phi(:,level:size(AQN.Y,2)-AQN.j-1)'*invMC(:,level+1:size(AQN.Y,2)-AQN.j+1));
  N = N(:,1);
  N = [zeros(level,1); N];

  % Run unit test for null space computation
  AQN.runUnitTest(4,invMC,rhs,level,phi,phi_comb,N);

  % Set coefficients of quadratic formula
  a = N'*invMC*N;
  b = 2*N'*invMC*phi_comb;
  c = phi_comb'*invMC*phi_comb - rhs(level,level);

  % Compute radicand
  radicand = b^2 - 4*a*c;

  % Check for real solution
  if radicand < 0

    % Print warning
    warning('AggQN: Imaginary root!  Radicand = %e',radicand);

    % Compute approximate lambda
    lambda = -b/(2*a);

  else

    % Compute lambda
    lambda = min(roots([a b c]));

  end

  % Compute "bottom" part component of Q*A
  AQN.MA_j(level+1:size(AQN.Y,2)-AQN.j+1,level) = lambda*N(level+1:size(AQN.Y,2)-AQN.j+1) + phi_comb(level+1:size(AQN.Y,2)-AQN.j+1) - AQN.S(:,level+AQN.j:end)'*(AQN.b_j(level)*AQN.y_j + AQN.Y(:,level+AQN.j-1));

end

% Compute invM*MA
AQN.A_j = invM*AQN.C'*AQN.MA_j;

% Run overall unit tests
AQN.runUnitTest(5);

end
