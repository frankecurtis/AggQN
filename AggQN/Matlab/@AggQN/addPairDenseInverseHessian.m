% Copyright (C) 2021 Albert Berahas, Frank E. Curtis, Baoyu Zhou
%
% All Rights Reserved.
%
% Authors: Albert Berahas, Frank E. Curtis, Baoyu Zhou
%
% Method definition for AggQN class

% Add pair, storage mode 'W'
function msg = addPairDenseInverseHessian(AQN,s,y)

% Compute inner product
r = 1/(s'*y);

% Compute intermediate vector
Wy = AQN.W*y;

% Compute inner product
yWy = y'*Wy;

% Update inverse Hessian (W)
AQN.W = AQN.W - r*Wy*s' - r*s*Wy' + r*(1 + r*yWy)*(s*s');

% Set message
msg = 'Add';

% Print message
if AQN.verbosity >= 1
  fprintf('AggQN: Pair added into inverse Hessian approximation.\n');
end

end
