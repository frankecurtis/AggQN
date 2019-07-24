% Copyright (C) 2019 Albert Berahas, Frank E. Curtis, Baoyu Zhou
%
% All Rights Reserved.
%
% Authors: Albert Berahas, Frank E. Curtis, Baoyu Zhou
%
% Method definition for AggQN class

% Add pair, storage mode 'H'
function msg = addPairH(AQN,s,y)

% Compute inner product
r = 1/(s'*y);

% Compute intermediate vector
Hs = AQN.H*s;

% Compute inner product
sHs = s'*Hs;

% Update Hessian (H)
AQN.H = AQN.H - (Hs*Hs')/sHs + r*(y*y');

% Set message
msg = 'Add';

% Print message
if AQN.verbosity >= 1
  fprintf('AggQN: Pair added into Hessian approximation.\n');
end

end