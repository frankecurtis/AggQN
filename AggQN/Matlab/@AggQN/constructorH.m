% Copyright (C) 2019 Albert Berahas, Frank E. Curtis, Baoyu Zhou
%
% All Rights Reserved.
%
% Authors: Albert Berahas, Frank E. Curtis, Baoyu Zhou
%
% Method definition for AggQN class

% Constructor for storage mode 'H'
function constructorH(AQN,H)

% Second input should be initial matrix
if ~isnumeric(H)
  msg = sprintf('AggQN: For storage mode ''H''');
  msg = [msg newline '       ... second input to constructor must be real symmetric matrix'];
  msg = [msg newline '       ... but given input is non-numeric'];
  msg = [msg newline AQN.usageMessage];
  error(msg);
end
if length(size(H)) ~= 2
  msg = sprintf('AggQN: For storage mode ''H''');
  msg = [msg newline '       ... second input to constructor must be real symmetric matrix'];
  msg = [msg newline '       ... but given input is not a two-dimensional matrix'];
  msg = [msg newline AQN.usageMessage];
  error(msg);
end
if size(H,1) ~= size(H,2)
  msg = sprintf('AggQN: For storage mode ''H''');
  msg = [msg newline '       ... second input to constructor must be real symmetric matrix'];
  msg = [msg newline '       ... but given input is not square'];
  msg = [msg newline AQN.usageMessage];
  error(msg);
end
if ~issymmetric(H)
  msg = sprintf('AggQN: For storage mode ''H''');
  msg = [msg newline '       ... second input to constructor must be real symmetric matrix'];
  msg = [msg newline '       ... but given input is not symmetric'];
  msg = [msg newline AQN.usageMessage];
  error(msg);
end

% Set initial matrix
AQN.H = H;
AQN.n = size(H,1);

% Reset values for full memory modes
AQN.adaptive  = false;
AQN.aggregate = false;
AQN.m         = inf;

end