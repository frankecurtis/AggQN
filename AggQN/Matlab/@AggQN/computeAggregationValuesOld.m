% Copyright (C) 2019 Albert Berahas, Frank E. Curtis, Baoyu Zhou
%
% All Rights Reserved.
%
% Authors: Albert Berahas, Frank E. Curtis, Baoyu Zhou
%
% Method definition for AggQN class

% Compute aggregation values 'A' and 'b'
function computeAggregationValuesOld(AQN)

% Construct H_{1:j-1} as AQN_H.hessian
H_matr = AQN.initHv(eye(AQN.n));
AQN_H  = AggQN('H',H_matr);
AQN_H.setVerbosity(0);
for i = 1:AQN.j-1
  AQN_H.addPair(AQN.S(:,i),AQN.Y(:,i));
end

% Construct HS_j
AQN.HS_j = AQN_H.hessian*AQN.S(:,AQN.j:end);

% Construct SHS_j
SHS_j = AQN.S(:,AQN.j:end)'*AQN.HS_j;

% Compute factorization of Q'*H*Q
[AQN.L_SHS,~] = AQN.choleskyPerturb(AQN.S(:,AQN.j:end)'*AQN.HS_j);

% Compute Q^{-1}
invM = AQN.L_SHS\eye(size(AQN.Y,2)-AQN.j+1);
invM = AQN.L_SHS'\invM;

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
AQN.runUnitTestOld(1,invM);

%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Compute "bottom" of Q*A %
%%%%%%%%%%%%%%%%%%%%%%%%%%%

% Compute Omega
Omega = AQN.S(:,AQN.j:end)'*AQN.y_j*AQN.b_j' + AQN.S(:,AQN.j:end)'*AQN.Y(:,AQN.j:end-1) - triu(AQN.SY(AQN.j:end,AQN.j:end-1),0);

% Compute omega
omega = AQN.b_j/sqrt(AQN.rho_j);

% Compute right-hand side
rhs = omega * omega' + Omega' * invM * Omega;

%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%
length = (size(AQN.Y,2)-AQN.j)*(size(AQN.Y,2)-AQN.j+1);
startPoint = zeros(length,1);
A_j_temp = AQN.mapToMatrix(startPoint,size(AQN.Y,2)-AQN.j);
vec_temp_1 = [];
for i = 1:size(A_j_temp,2)
    vec_temp_1 = [vec_temp_1 ; SHS_j(1:i,:)*A_j_temp(:,i) + AQN.b_j(i)*AQN.S(:,1:i)'*AQN.y_j];
end

func_temp_2 = A_j_temp' * SHS_j * A_j_temp + Omega'*A_j_temp + A_j_temp'*Omega - omega*omega';
vec_temp_2 = AQN.mapToVector(func_temp_2);

vec_temp = [vec_temp_1 ; vec_temp_2];

% Evaluate start norm of vec_temp
AQN.startValue = norm(vec_temp,inf);
%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%

% Solve for the first quadratic equation (eq. (3.21c)), solved directly
% without computing a_{m-1,2}^* + lambda_{m-1}
AQN.MA_j(size(AQN.Y,2)-AQN.j+1,size(AQN.Y,2)-AQN.j) = sqrt(rhs(size(AQN.Y,2)-AQN.j,size(AQN.Y,2)-AQN.j)/invM(size(AQN.Y,2)-AQN.j+1,size(AQN.Y,2)-AQN.j+1)) - AQN.S(:,end)'*(AQN.b_j(size(AQN.Y,2)-AQN.j)*AQN.y_j + AQN.Y(:,size(AQN.Y,2)-1));

% Call unit test for first quadratic
AQN.runUnitTestOld(2,invM,vec_temp,rhs);

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
    beta_sub = (phi(:,good_ind)'*invM*phi(:,good_ind))\rhs(good_ind+1,level);
    beta     = zeros(size(AQN.Y,2)-AQN.j-level,1);
    for i = 1:size(good_ind,1)
      beta(good_ind(i)+1-level) = beta_sub(i);
    end
  else
    good_ind = good_ind(2:end);
    beta_sub = (phi(:,good_ind)'*invM*phi(:,good_ind))\rhs(good_ind+1,level);
    beta     = zeros(size(AQN.Y,2)-AQN.j-level,1);
    for i = 1:size(good_ind,1)
      beta(good_ind(i)+1-level) = beta_sub(i);
    end
  end

  % Compute combination of phi
  phi_comb = phi(:,level:size(AQN.Y,2)-AQN.j-1)*beta;

  % Call unit test for phi combination
  AQN.runUnitTestOld(3,invM,vec_temp,rhs,level,phi,phi_comb);
  
  % Compute null space for affine equations
  N = null(phi(:,level:size(AQN.Y,2)-AQN.j-1)'*invM(:,level+1:size(AQN.Y,2)-AQN.j+1));
  N = N(:,1);
  N = [zeros(level,1); N];

  % Run unit test for null space computation
  AQN.runUnitTestOld(4,invM,vec_temp,rhs,level,phi,phi_comb,N);

  % Set coefficients of quadratic formula
  a = N'*invM*N;
  b = 2*N'*invM*phi_comb;
  c = phi_comb'*invM*phi_comb - rhs(level,level);

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
AQN.A_j = invM*AQN.MA_j;

%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%
vec_temp_1 = [];
for i = 1:size(A_j_temp,2)
    vec_temp_1 = [vec_temp_1 ; SHS_j(1:i,:)*AQN.A_j(:,i) + AQN.b_j(i)*AQN.S(:,1:i)'*AQN.y_j];
end

func_temp_2 = AQN.A_j' * SHS_j * AQN.A_j + Omega'*AQN.A_j + AQN.A_j'*Omega - omega*omega';
vec_temp_2 = AQN.mapToVector(func_temp_2);

vec_temp = [vec_temp_1 ; vec_temp_2];

AQN.checkFlag = AQN.checkAccuracy(vec_temp);

AQN.DSacc = norm(vec_temp,inf)/max(1,AQN.startValue);
%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%

% Run overall unit tests
AQN.runUnitTestOld(5,invM,vec_temp);

    

end
