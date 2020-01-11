% Copyright (C) 2019 Albert Berahas, Frank E. Curtis, Baoyu Zhou
%
% All Rights Reserved.
%
% Authors: Albert Berahas, Frank E. Curtis, Baoyu Zhou
%
% Method definition for AggQN class

% Compute aggregation values 'A' and 'b'
function computeAggregationValuesNewton(AQN)

% Set length
length = size(AQN.Y,2)-AQN.j;

% Set starting point
vecA_j = zeros(length*(length+1),1);

% Set starting point as matrix
A_j = reshape(vecA_j,[length+1,length]);

% Evaluate error and set temp vectors
[errorMaxAbs,errorFroSqr,v1,v2] = AQN.computeAggregationValuesError(A_j);

% Initialize iteration number
iteration_counter = 0;

% Initialize regularization
regularization = AQN.reg_newton;

% Print header
if AQN.verbosity >= 2
  fprintf('%6d %13s %13s %13s %13s %13s %13s %13s %13s %13s %13s',...
          'iter.','errorMaxAbs','errorFroSqr',...
          '|d Newton|','dd Newton','error Newton','alpha Newton',...
          '|d steep|','dd steep','error steep','alpha steep');
end

% Iteration loop
while 1

  % Print iteration
  if AQN.verbosity >= 2
    fprintf('%6d %e %e',iteration_counter,errorMaxAbs,errorFroSqr);
  end

  % Check error
  if errorMaxAbs <= AQN.acc_tol_newton
    break;
  elseif iteration_counter > AQN.max_iter_newton
    break;
  end
  
  % Compute Jacobian
  Jacobian = [];
  for i = 1:(length)
    Jacobian = [Jacobian ; [zeros(i,(length+1)*(i-1)) AQN.SHS_j(1:i,:) zeros(i,(length+1)*(length-i))]];
  end
  for i = 1:(length)
    for j = 1:i
      row1 = [zeros(1,(length+1)*(j-1)) (A_j(:,i)'*AQN.SHS_j + AQN.Omega_j(:,i)') zeros(1,(length+1)*((length)-j))];
      row2 = [zeros(1,(length+1)*(i-1)) (A_j(:,j)'*AQN.SHS_j + AQN.Omega_j(:,j)') zeros(1,(length+1)*((length)-i))];
      Jacobian = [Jacobian ; row1 + row2];
    end
  end
  
  % Compute Newton and steepest descent steps
  errorFroSqr_trial = inf*ones(2,1);
  d = zeros(size(Jacobian,1),2);
  alpha = zeros(2,1);
  for step_type = 1:1
  
    % Compute direction
    if step_type == 1
      d(:,step_type) = -(Jacobian'*Jacobian + regularization*eye(size(Jacobian,2)))\(Jacobian'*[v1;v2]);
    else
      d(:,step_type) = -Jacobian'*[v1;v2];
    end
    
    % Print step
    if AQN.verbosity >= 2
      fprintf(' %e %+e',norm(d(:,step_type)),[v1;v2]'*Jacobian*d(:,step_type));
    end
  
    % Initialize stepsize
    alpha(step_type) = 1;
  
    % Line search
    while 1
    
      % Set trial point
      vecA_j_trial = vecA_j + alpha(step_type)*d(:,step_type);
      
      % Set trial point as matrix
      A_j_trial = reshape(vecA_j_trial,[length+1,length]);
    
      % Evaluate trial objective
      [~,errorFroSqr_trial(step_type),~,~] = AQN.computeAggregationValuesError(A_j_trial);
    
      % Check for sufficient decrease
      if errorFroSqr_trial(step_type) <= errorFroSqr + alpha*AQN.c1_newton*([v1;v2]'*Jacobian*d(:,step_type))
        break;
      elseif alpha(step_type) <= AQN.alpha_min_newton
        break;
      else
        alpha(step_type) = alpha(step_type)/2;
      end
    
    end
    
    % Print trial information
    if AQN.verbosity >= 2
      fprintf(' %e %e',alpha(step_type),errorFroSqr_trial(step_type));
    end
  
  end
  
  % Print new line
  if AQN.verbosity >= 2
    fprintf('\n');
  end
  
  % Choose step
  if errorFroSqr_trial(1) <= errorFroSqr_trial(2)
    vecA_j = vecA_j + alpha(1)*d(:,1);
  else
    vecA_j = vecA_j + alpha(2)*d(:,2);
  end
  A_j = reshape(vecA_j,[length+1,length]);
  
  % Evaluate error and set temp vectors
  [errorMaxAbs,errorFroSqr,v1,v2] = AQN.computeAggregationValuesError(A_j);
  
  % Increment iteration counter
  iteration_counter = iteration_counter + 1;
  
end

% Set result
AQN.A_j = A_j;

end