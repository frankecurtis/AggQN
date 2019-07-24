% Copyright (C) 2019 Albert Berahas, Frank E. Curtis, Baoyu Zhou
%
% All Rights Reserved.
%
% Authors: Albert Berahas, Frank E. Curtis, Baoyu Zhou
%
% Method definition for AggQN class

function addDataSY(AQN,s,y)

% Add to S*Y matrix
if size(AQN.S,2) >= 1
  AQN.SY = [AQN.SY (y'*AQN.S)'; s'*AQN.Y s'*y];
else
  AQN.SY = s'*y;
end

% Add to S and Y
AQN.S = [AQN.S s];
AQN.Y = [AQN.Y y];

% Add to rho
AQN.rho = [AQN.rho; 1/(AQN.SY(end,end))];

% Factorizations needed when aggregating
if AQN.aggregate
  
  %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
  % UPDATE CHOLESKY FOR S'*S (REVERSE ORDER) %
  %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
  
  % Reverse order
  S_reverse = zeros(size(AQN.S));
  for i = 1:size(AQN.S,2)
    S_reverse(:,i) = AQN.S(:,end+1-i);
  end
  
  % Compute factorization
  [AQN.L_SS,perturbation] = AQN.choleskyPerturb(S_reverse'*S_reverse);
  
  % Check for error
  if AQN.verbosity >= 1
    fprintf('AggQN: Cholesky factorization of S_reverse''*S_reverse performed, perturbation = %e\n',perturbation);
  end
    
end

end