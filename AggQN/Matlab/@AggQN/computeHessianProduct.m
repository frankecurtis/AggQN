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
if strcmp(AQN.storage_mode,'SY') == 1
  
  % TO DO
  error('AggQN: Hessian-vector product not yet implemented for storage mode ''SY''.');
  
elseif strcmp(AQN.storage_mode,'W') == 1
  
  % Solve system
  Hv = AQN.W\v;
  
  % Print message
  if AQN.verbosity >= 1
    fprintf('AggQN: Computed Hessian-vector product.  Inefficient since storage mode is ''W''.\n');
  end
  
else % strcmp(AQN.storage_mode,'H') == 1
  
  % Compute product
  Hv = AQN.H*v;
  
  % Print message
  if AQN.verbosity >= 1
    fprintf('AggQN: Computed Hessian-vector product.\n');
  end
  
end

end