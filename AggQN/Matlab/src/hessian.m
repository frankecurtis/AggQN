% Copyright (C) 2019 Albert Berahas, Frank E. Curtis, Baoyu Zhou
%
% All Rights Reserved.
%
% Authors: Albert Berahas, Frank E. Curtis, Baoyu Zhou
%
% Method definition for AggQN class

% Hessian
function H = hessian(AQN)

% Check storage mode
if strcmp(AQN.storage_mode,'SY') == 1

  % Set return value
  H = inv(AQN.inverseHessian);
  
  % Print message
  warning('AggQN: Constructing Hessian for storage mode ''W'' is inefficient and inaccurate!');
  
elseif strcmp(AQN.storage_mode,'W') == 1
  
  % Set return value
  H = inv(AQN.W);
  
  % Print message
  warning('AggQN: Constructing Hessian for storage mode ''W'' is inefficient and inaccurate!');
  
else % strcmp(AQN.storage_mode,'H') == 1
  
  % Set return value
  H = AQN.H;
  
end

end
