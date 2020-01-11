% Copyright (C) 2019 Albert Berahas, Frank E. Curtis, Baoyu Zhou
%
% All Rights Reserved.
%
% Authors: Albert Berahas, Frank E. Curtis, Baoyu Zhou
%
% Method definition for AggQN class

% Sets verbosity
function setAccuracyTolerance(AQN,value)

% Check input type
if length(value) ~= 1 || ~isnumeric(value) || value < 0.0
  error('AggQN: Invalid input to setAccuracyTolerance.  Input must be nonnegative scalar.');
end

% Print message
if AQN.acc_tol ~= value && AQN.verbosity >= 1
  fprintf('AggQN: Accuracy tolerance set to %d\n',value);
end

% If all is OK, then set option
AQN.acc_tol = value;

end