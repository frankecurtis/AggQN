% Copyright (C) 2021 Albert Berahas, Frank E. Curtis, Baoyu Zhou
%
% All Rights Reserved.
%
% Authors: Albert Berahas, Frank E. Curtis, Baoyu Zhou
%
% Method definition for AggQN class

% Hessian
function H = hessian(AQN)

% Check storage mode
if strcmp(AQN.storage_mode,'limitedMemory') == 1

  % Set return value (using compact form)
  if size(AQN.S,2) >= 1
    H = AQN.initHv(eye(AQN.n)) - [AQN.HS AQN.Y]*([AQN.SHS AQN.L; AQN.L' -AQN.D]\[AQN.HS'; AQN.Y']);
  else
    H = AQN.initHv(eye(AQN.n));
  end

  % Print message
  warning('AggQN: Constructing Hessian for storage mode ''limitedMemory'' is inefficient and inaccurate!');

elseif strcmp(AQN.storage_mode,'denseInverseHessian') == 1

  % Set return value
  H = inv(AQN.W);

  % Print message
  warning('AggQN: Constructing Hessian for storage mode ''denseInverseHessian'' is inefficient and inaccurate!');

else % strcmp(AQN.storage_mode,'denseHessian') == 1

  % Set return value
  H = AQN.H;

end

end
