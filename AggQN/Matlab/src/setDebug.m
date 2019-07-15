% Copyright (C) 2019 Albert Berahas, Frank E. Curtis, Baoyu Zhou
%
% All Rights Reserved.
%
% Authors: Albert Berahas, Frank E. Curtis, Baoyu Zhou
%
% Method definition for AggQN class

% Sets debug
function setDebug(AQN,value)

% Check input type
if length(value) ~= 1 || ~islogical(value)
  error('AggQN: Invalid input to setDebug.  Input must be logical scalar (true or false).');
end

% Print message
if AQN.debug ~= value && AQN.verbosity >= 1
  fprintf('AggQN: Debug set to %d\n',value);
end

% If all is OK, then set option
AQN.debug = value;

% If all is OK, then set verbosity to at least 1
if AQN.verbosity < 1
  AQN.setVerbosity(1);
end

end