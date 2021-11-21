% Copyright (C) 2021 Albert Berahas, Frank E. Curtis, Baoyu Zhou
%
% All Rights Reserved.
%
% Authors: Albert Berahas, Frank E. Curtis, Baoyu Zhou
%
% Method definition for AggQN class

% Compute inverse-Hessian-vector product (W*v)
function Wv = computeInverseHessianProduct(AQN,v)

% Check option
if strcmp(AQN.storage_mode,'limitedMemory') == 1

  % Check if pairs exist
  if size(AQN.S,2) == 0
    Wv = AQN.initWv(v);
  else

    % Do two-loop recursion
    alpha = zeros(size(AQN.S,2),1);
    for i = size(AQN.S,2):-1:1
      alpha(i) = AQN.rho(i)*AQN.S(:,i)'*v;
      v        = v - alpha(i)*AQN.Y(:,i);
    end
    Wv = AQN.initWv(v);
    for i = 1:size(AQN.S,2)
      beta = AQN.rho(i)*AQN.Y(:,i)'*Wv;
      Wv   = Wv + AQN.S(:,i)*(alpha(i)-beta);
    end

  end

  % Print message
  if AQN.verbosity >= 1
    fprintf('AggQN: Computed inverse-Hessian-vector product by two-loop recursion.\n');
  end

elseif strcmp(AQN.storage_mode,'denseInverseHessian') == 1

  % Compute product
  Wv = AQN.W*v;

  % Print message
  if AQN.verbosity >= 1
    fprintf('AggQN: Computed inverse-Hessian-vector product.\n');
  end

else % strcmp(AQN.storage_mode,'denseHessian') == 1

  % Solve system
  Wv = AQN.H\v;

  % Print message
  if AQN.verbosity >= 1
    fprintf('AggQN: Computed inverse-Hessian-vector product.  Inefficient since storage mode is ''denseHessian''.\n');
  end

end

end
