% Copyright (C) 2021 Albert Berahas, Frank E. Curtis, Baoyu Zhou
%
% All Rights Reserved.
%
% Authors: Albert Berahas, Frank E. Curtis, Baoyu Zhou
%
% Method definition for AggQN class

% Compute Hessian-vector product (H*v)
function Hv = computeHessianProduct(AQN,v)

% Check option
if strcmp(AQN.storage_mode,'limitedMemory') == 1

  % Check if pairs exist
  if size(AQN.S,2) >= 1

    % Compute compact form product
    temp = [AQN.SHS AQN.L; AQN.L' -AQN.D]\([AQN.HS'; AQN.Y']*v);

    % Compute product
    Hv = AQN.initHv(v) - [AQN.HS AQN.Y]*temp;

  else

    % Compute product
    Hv = AQN.initHv(v);

  end

  % Print message
  if AQN.verbosity >= 1
    fprintf('AggQN: Computed Hessian-vector product by compact form.\n');
  end

elseif strcmp(AQN.storage_mode,'denseInverseHessian') == 1

  % Solve system
  Hv = AQN.W\v;

  % Print message
  if AQN.verbosity >= 1
    fprintf('AggQN: Computed Hessian-vector product.  Inefficient since storage mode is ''denseInverseHessian''.\n');
  end

else % strcmp(AQN.storage_mode,'denseHessian') == 1

  % Compute product
  Hv = AQN.H*v;

  % Print message
  if AQN.verbosity >= 1
    fprintf('AggQN: Computed Hessian-vector product.\n');
  end

end

end
