% Copyright (C) 2019 Albert Berahas, Frank E. Curtis, Baoyu Zhou
%
% All Rights Reserved.
%
% Authors: Albert Berahas, Frank E. Curtis, Baoyu Zhou
%
% Method definition for AggQN class

% Constructor for storage mode 'W'
function constructorW(AQN,W)

% Second input should be initial matrix
if ~isnumeric(W)
  msg = sprintf('AggQN: For storage mode ''W''');
  msg = [msg newline '       ... second input to constructor must be real symmetric matrix'];
  msg = [msg newline '       ... but given input is non-numeric'];
  msg = [msg newline AQN.usageMessage];
  error(msg);
end
if length(size(W)) ~= 2
  msg = sprintf('AggQN: For storage mode ''W''');
  msg = [msg newline '       ... second input to constructor must be real symmetric matrix'];
  msg = [msg newline '       ... but given input is not a two-dimensional matrix'];
  msg = [msg newline AQN.usageMessage];
  error(msg);
end
if size(W,1) ~= size(W,2)
  msg = sprintf('AggQN: For storage mode ''W''');
  msg = [msg newline '       ... second input to constructor must be real symmetric matrix'];
  msg = [msg newline '       ... but given input is not square'];
  msg = [msg newline AQN.usageMessage];
  error(msg);
end
if ~issymmetric(W)
  msg = sprintf('AggQN: For storage mode ''W''');
  msg = [msg newline '       ... second input to constructor must be real symmetric matrix'];
  msg = [msg newline '       ... but given input is not symmetric'];
  msg = [msg newline AQN.usageMessage];
  error(msg);
end

% Set initial matrix
AQN.W = W;
AQN.n = size(W,1);

% Reset values for full memory modes
AQN.adaptive  = false;
AQN.aggregate = false;
AQN.m         = inf;

end