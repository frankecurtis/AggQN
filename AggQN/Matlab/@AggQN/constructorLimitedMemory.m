% Copyright (C) 2021 Albert Berahas, Frank E. Curtis, Baoyu Zhou
%
% All Rights Reserved.
%
% Authors: Albert Berahas, Frank E. Curtis, Baoyu Zhou
%
% Method definition for AggQN class

% Constructor for storage mode 'limitedMemory'
function constructorLimitedMemory(AQN,initWv,initHv,n,m)

% Second and third inputs should be function handles
if ~isa(initWv,'function_handle') || ~isa(initHv,'function_handle')
  msg = 'AggQN: For storage mode ''limitedMemory''';
  msg = [msg newline '       ... 2nd input to constructor must be function handle for computing initial inverse-Hessian-vector product,'];
  msg = [msg newline '       ... 3rd input to constructor must be function handle for computing initial Hessian-vector product,'];
  msg = [msg newline '       ... 4th input to constructor must be vector length, and'];
  msg = [msg newline '       ... 5th input to constructor must be history length.'];
  msg = [msg newline AQN.usageMessage];
  error(msg);
end

% Set initial matrix
AQN.initWv = initWv;
AQN.initHv = initHv;

% Fourth input must be vector length, a positive integer
if ~isnumeric(n)
  msg = 'AggQN: For storage mode ''limitedMemory''';
  msg = [msg newline '       ... 4th input to constructor must be positive integer'];
  msg = [msg newline '       ... but given input is non-numeric'];
  msg = [msg newline AQN.usageMessage];
  error(msg);
end
if length(n) ~= 1
  msg = 'AggQN: For storage mode ''limitedMemory''';
  msg = [msg newline '       ... 4th input to constructor must be positive integer'];
  msg = [msg newline '       ... but given input is not a scalar'];
  msg = [msg newline AQN.usageMessage];
  error(msg);
end
if n ~= abs(n) || n < 1 || n == inf
  msg = 'AggQN: For storage mode ''limitedMemory''';
  msg = [msg newline '       ... 4th input to constructor must be positive integer'];
  msg = [msg newline AQN.usageMessage];
  error(msg);
end

% Set vector length
AQN.n = n;

% Fifth input must be history, a positive integer
if ~isnumeric(m)
  msg = 'AggQN: For storage mode ''limitedMemory''';
  msg = [msg newline '       ... 5th input to constructor must be positive integer or inf'];
  msg = [msg newline '       ... but given input is non-numeric'];
  msg = [msg newline AQN.usageMessage];
  error(msg);
end
if length(m) ~= 1
  msg = 'AggQN: For storage mode ''limitedMemory''';
  msg = [msg newline '       ... 5th input to constructor must be positive integer or inf'];
  msg = [msg newline '       ... but given input is not a scalar'];
  msg = [msg newline AQN.usageMessage];
  error(msg);
end
if m ~= abs(m) || m < 1
  msg = 'AggQN: For storage mode ''limitedMemory''';
  msg = [msg newline '       ... 5th input to constructor must be positive integer or inf'];
  msg = [msg newline AQN.usageMessage];
  error(msg);
end

% Set history length
AQN.m = m;

end
