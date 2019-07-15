% Copyright (C) 2019 Albert Berahas, Frank E. Curtis, Baoyu Zhou
%
% All Rights Reserved.
%
% Authors: Albert Berahas, Frank E. Curtis, Baoyu Zhou
%
% Method definition for AggQN class

% Constructor for storage mode 'SY'
function constructorSY(AQN,initWv,initHv,n,m)

% Second and third inputs should be function handles
if ~isa(initWv,'function_handle') || ~isa(initHv,'function_handle')
  msg = sprintf('AggQN: For storage mode ''SY''');
  msg = [msg newline '       ... second input to constructor must be function handle'];
  msg = [msg newline '       ... for computing initial inverse-Hessian-vector product'];
  msg = [msg newline '       ... and third input to constructor must be function handle'];
  msg = [msg newline '       ... for computing initial Hessian-vector product'];
  msg = [msg newline AQN.usageMessage];
  error(msg);
end

% Set initial matrix
AQN.initWv = initWv;
AQN.initHv = initHv;

% Fourth input must be vector length, a positive integer
if ~isnumeric(n)
  msg = 'AggQN: Fourth input to constructor must be positive integer';
  msg = [msg newline '       ... but given input is non-numeric'];
  msg = [msg newline AQN.usageMessage];
  error(msg);
end
if length(n) ~= 1
  msg = 'AggQN: Fourth input to constructor must be positive integer';
  msg = [msg newline '       ... but given input is not a scalar'];
  msg = [msg newline AQN.usageMessage];
  error(msg);
end
if n ~= abs(n) || n < 1 || n == inf
  msg = 'AggQN: Fourth input to constructor must be positive integer';
  msg = [msg newline AQN.usageMessage];
  error(msg);
end

% Set vector length
AQN.n = n;

% Fifth input must be history, a positive integer
if ~isnumeric(m)
  msg = 'AggQN: Fifth input to constructor must be positive integer or inf';
  msg = [msg newline '       ... but given input is non-numeric'];
  msg = [msg newline AQN.usageMessage];
  error(msg);
end
if length(m) ~= 1
  msg = 'AggQN: Fifth input to constructor must be positive integer or inf';
  msg = [msg newline '       ... but given input is not a scalar'];
  msg = [msg newline AQN.usageMessage];
  error(msg);
end
if m ~= abs(m) || m < 1
  msg = 'AggQN: Fifth input to constructor must be positive integer or inf';
  msg = [msg newline AQN.usageMessage];
  error(msg);
end

% Set history
AQN.m = m;

end