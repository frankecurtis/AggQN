% Copyright (C) 2021 Albert Berahas, Frank E. Curtis, Baoyu Zhou
%
% All Rights Reserved.
%
% Authors: Albert Berahas, Frank E. Curtis, Baoyu Zhou
%
% Method definition for AggQN class

% Compute "bottom" part of A with exact method
function computeAggregationValuesExact(AQN)

% Solve for the first quadratic equation (eq. (3.21c)), solved directly
% without computing a_{m-1,2}^* + lambda_{m-1}
AQN.MA_j(size(AQN.Y,2)-AQN.j+1,size(AQN.Y,2)-AQN.j) = sqrt(AQN.rhs_j(size(AQN.Y,2)-AQN.j,size(AQN.Y,2)-AQN.j)/AQN.invM_j(size(AQN.Y,2)-AQN.j+1,size(AQN.Y,2)-AQN.j+1)) - AQN.S(:,end)'*(AQN.b_j(size(AQN.Y,2)-AQN.j)*AQN.y_j + AQN.Y(:,size(AQN.Y,2)-1));

% Call unit test for first quadratic
AQN.runUnitTest(2);

% Initialize phi
phi = zeros(size(AQN.Y,2)-AQN.j+1,size(AQN.Y,2)-AQN.j-1);

% Initialize values for "cutting off" beta
good_ind = [];

% Loop through remaining levels in reverse order
for level = size(AQN.Y,2)-AQN.j-1:-1:1

  % Compute phi
  phi(level+2:size(AQN.Y,2)-AQN.j+1,level) = AQN.MA_j(level+2:size(AQN.Y,2)-AQN.j+1,level+1) + AQN.S(:,level+AQN.j+1:end)'*(AQN.b_j(level+1)*AQN.y_j + AQN.Y(:,level+AQN.j));

  % Determine if new column of phi is in span of previous ones
  tau        = (phi(:,good_ind)'*phi(:,good_ind))\(phi(:,good_ind)'*phi(:,level));
  projection = phi(:,good_ind)*tau;
  error      = phi(:,level) - projection;

  % Compute beta
  if norm(error) > AQN.proj_tol_beta * norm(projection)
    good_ind = [level; good_ind];
    beta_sub = (phi(:,good_ind)'*AQN.invM_j*phi(:,good_ind))\AQN.rhs_j(good_ind+1,level);
    beta     = zeros(size(AQN.Y,2)-AQN.j-level,1);
    for i = 1:size(good_ind,1)
      beta(good_ind(i)+1-level) = beta_sub(i);
    end
  else
    beta_sub = (phi(:,good_ind)'*AQN.invM_j*phi(:,good_ind))\AQN.rhs_j(good_ind+1,level);
    beta     = zeros(size(AQN.Y,2)-AQN.j-level,1);
    for i = 1:size(good_ind,1)
      beta(good_ind(i)+1-level) = beta_sub(i);
    end
  end

  % Compute combination of phi
  phi_comb = phi(:,level:size(AQN.Y,2)-AQN.j-1)*beta;

  % Call unit test for phi combination
  AQN.runUnitTest(3,level,phi,phi_comb);

  % Compute null space for affine equations
  N = null(phi(:,level:size(AQN.Y,2)-AQN.j-1)'*AQN.invM_j(:,level+1:size(AQN.Y,2)-AQN.j+1));
  N = N(:,1);
  N = [zeros(level,1); N];

  % Run unit test for null space computation
  AQN.runUnitTest(4,level,phi,phi_comb,N);

  % Set coefficients of quadratic formula
  a = N'*AQN.invM_j*N;
  b = 2*N'*AQN.invM_j*phi_comb;
  c = phi_comb'*AQN.invM_j*phi_comb - AQN.rhs_j(level,level);

  % Compute radicand
  radicand = b^2 - 4*a*c;

  % Check for real solution
  if radicand < 0

    % Print warning
    warning('AggQN: Imaginary root!  Radicand = %e',radicand);

    % Compute approximate lambda
    lambda = -b/(2*a);

  elseif isnan(radicand) || isinf(radicand)

    % Abort!
    break;

  else

    % Compute lambda
    lambda = min(roots([a b c]));

  end

  % Compute "bottom" part component of Q*A
  AQN.MA_j(level+1:size(AQN.Y,2)-AQN.j+1,level) = lambda*N(level+1:size(AQN.Y,2)-AQN.j+1) + phi_comb(level+1:size(AQN.Y,2)-AQN.j+1) - AQN.S(:,level+AQN.j:end)'*(AQN.b_j(level)*AQN.y_j + AQN.Y(:,level+AQN.j-1));

end

% Compute invM_j*MA
AQN.A_j = AQN.invM_j*AQN.MA_j;

end
