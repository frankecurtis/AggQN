% Copyright (C) 2021 Albert Berahas, Frank E. Curtis, Baoyu Zhou
%
% All Rights Reserved.
%
% Authors: Albert Berahas, Frank E. Curtis, Baoyu Zhou
%
% Method definition for AggQN class

% Inverse Hessian
function W = inverseHessian(AQN)

% Check storage mode
if strcmp(AQN.storage_mode,'limitedMemory') == 1

  % Check size
  if size(AQN.S,2) >= 1

    % Construct intermediate matrices
    RinvSt   = AQN.R\(AQN.S');
    YRinvSt  = AQN.Y*RinvSt;
    WYRinvSt = AQN.initWv(YRinvSt);

    % Construct inverse Hessian
    W = AQN.initWv(eye(AQN.n)) - WYRinvSt - WYRinvSt' + RinvSt'*AQN.D*RinvSt + YRinvSt'*WYRinvSt;

  else

    % Construct inverse Hessian
    W = AQN.initWv(eye(AQN.n));

  end

  % Print message
  warning('AggQN: Constructing inverse Hessian for storage mode ''limitedMemory'' is inefficient and inaccurate!');

elseif strcmp(AQN.storage_mode,'denseInverseHessian') == 1

  % Set return value
  W = AQN.W;

else % strcmp(AQN.storage_mode,'denseHessian') == 1

  % Set return value
  W = inv(AQN.H);

  % Print message
  warning('AggQN: Constructing inverse Hessian for storage mode ''denseHessian'' is inefficient and inaccurate!');

end

end
