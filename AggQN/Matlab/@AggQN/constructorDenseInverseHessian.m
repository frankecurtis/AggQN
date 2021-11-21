% Copyright (C) 2021 Albert Berahas, Frank E. Curtis, Baoyu Zhou
%
% All Rights Reserved.
%
% Authors: Albert Berahas, Frank E. Curtis, Baoyu Zhou
%
% Method definition for AggQN class

% Constructor for storage mode 'denseInverseHessian'
function constructorDenseInverseHessian(AQN,W)

% Second input should be initial matrix
if ~isnumeric(W)
  msg = 'AggQN: For storage mode ''denseInverseHessian''';
  msg = [msg newline '       ... 2nd input to constructor must be real symmetric matrix'];
  msg = [msg newline '       ... but given input is non-numeric'];
  msg = [msg newline AQN.usageMessage];
  error(msg);
end
if length(size(W)) ~= 2
  msg = 'AggQN: For storage mode ''denseInverseHessian''';
  msg = [msg newline '       ... 2nd input to constructor must be real symmetric matrix'];
  msg = [msg newline '       ... but given input is not a two-dimensional matrix'];
  msg = [msg newline AQN.usageMessage];
  error(msg);
end
if size(W,1) ~= size(W,2)
  msg = 'AggQN: For storage mode ''denseInverseHessian''';
  msg = [msg newline '       ... 2nd input to constructor must be real symmetric matrix'];
  msg = [msg newline '       ... but given input is not square'];
  msg = [msg newline AQN.usageMessage];
  error(msg);
end
if ~issymmetric(W)
  msg = 'AggQN: For storage mode ''denseInverseHessian''';
  msg = [msg newline '       ... 2nd input to constructor must be real symmetric matrix'];
  msg = [msg newline '       ... but given input is not symmetric'];
  msg = [msg newline AQN.usageMessage];
  error(msg);
end

% Set initial matrix
AQN.W = W;
AQN.n = size(W,1);

% Reset values for full memory modes
AQN.aggregate = false;
AQN.onlySwap  = false;
AQN.tryNewton = false;
AQN.m         = inf;

end
