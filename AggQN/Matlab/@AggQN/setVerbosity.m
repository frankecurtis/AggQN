% Copyright (C) 2021 Albert Berahas, Frank E. Curtis, Baoyu Zhou
%
% All Rights Reserved.
%
% Authors: Albert Berahas, Frank E. Curtis, Baoyu Zhou
%
% Method definition for AggQN class

% Sets verbosity
function setVerbosity(AQN,value)

% Check input type
if length(value) ~= 1 || ~isnumeric(value) || value ~= abs(value)
  error('AggQN: Invalid input to setVerbosity.  Input must be nonnegative integer.');
end

% Print message
if AQN.verbosity ~= value && value >= 1
  fprintf('AggQN: Verbosity set to %d\n',value);
end

% If all is OK, then set option
AQN.verbosity = value;

end
