% Copyright (C) 2019 Albert Berahas, Frank E. Curtis, Baoyu Zhou
%
% All Rights Reserved.
%
% Authors: Albert Berahas, Frank E. Curtis, Baoyu Zhou
%
% Method definition for AggQN class

% Constructor for storage mode 'denseHessian'
function constructorDenseHessian(AQN,H)

% Second input should be initial matrix
if ~isnumeric(H)
  msg = 'AggQN: For storage mode ''denseHessian''';
  msg = [msg newline '       ... 2nd input to constructor must be real symmetric matrix'];
  msg = [msg newline '       ... but given input is non-numeric'];
  msg = [msg newline AQN.usageMessage];
  error(msg);
end
if length(size(H)) ~= 2
  msg = 'AggQN: For storage mode ''denseHessian''';
  msg = [msg newline '       ... 2nd input to constructor must be real symmetric matrix'];
  msg = [msg newline '       ... but given input is not a two-dimensional matrix'];
  msg = [msg newline AQN.usageMessage];
  error(msg);
end
if size(H,1) ~= size(H,2)
  msg = 'AggQN: For storage mode ''denseHessian''';
  msg = [msg newline '       ... 2nd input to constructor must be real symmetric matrix'];
  msg = [msg newline '       ... but given input is not square'];
  msg = [msg newline AQN.usageMessage];
  error(msg);
end
if ~issymmetric(H)
  msg = 'AggQN: For storage mode ''denseHessian''';
  msg = [msg newline '       ... 2nd input to constructor must be real symmetric matrix'];
  msg = [msg newline '       ... but given input is not symmetric'];
  msg = [msg newline AQN.usageMessage];
  error(msg);
end

% Set initial matrix
AQN.H = H;
AQN.n = size(H,1);

% Reset values for full memory modes
AQN.aggregate = false;
AQN.onlySwap  = false;
AQN.tryNewton = false;
AQN.m         = inf;

end