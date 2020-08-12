% Copyright (C) 2019 Albert Berahas, Frank E. Curtis, Baoyu Zhou
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
  
  % TO DO
  error('AggQN: Hessian-vector product not yet implemented for storage mode ''limitedMemory''.');
  
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