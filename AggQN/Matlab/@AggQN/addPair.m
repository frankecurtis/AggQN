% Copyright (C) 2021 Albert Berahas, Frank E. Curtis, Baoyu Zhou
%
% All Rights Reserved.
%
% Authors: Albert Berahas, Frank E. Curtis, Baoyu Zhou
%
% Method definition for AggQN class

% Add pair
function msg = addPair(AQN,s,y)

% Check size
if length(size(s)) ~= 2     || ...
   size(s,1)       ~= AQN.n || ...
   size(s,2)       ~= 1     || ...
   length(size(y)) ~= 2     || ...
   size(y,1)       ~= AQN.n || ...
   size(y,2)       ~= 1
  error('AggQN: Invalid input to addPair(s,y).  Inputs s and y must be column vectors of length %d for this object.',AQN.n);
end

% Check storage mode
if strcmp(AQN.storage_mode,'limitedMemory') == 1

  % Add pair
  msg = AQN.addPairLimitedMemory(s,y);

elseif strcmp(AQN.storage_mode,'denseInverseHessian') == 1

  % Add pair
  msg = AQN.addPairDenseInverseHessian(s,y);

else % strcmp(AQN.storage_mode,'denseHessian') == 1

  % Add pair
  msg = AQN.addPairDenseHessian(s,y);

end

% Check verbosity
if AQN.verbosity >= 2, AQN.printData; end

end
