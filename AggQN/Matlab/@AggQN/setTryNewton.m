% Copyright (C) 2021 Albert Berahas, Frank E. Curtis, Baoyu Zhou
%
% All Rights Reserved.
%
% Authors: Albert Berahas, Frank E. Curtis, Baoyu Zhou
%
% Method definition for AggQN class

% Sets verbosity
function setTryNewton(AQN,value)

% Check input type
if length(value) ~= 1 || ~islogical(value)
  error('AggQN: Invalid input to setTryNewton.  Input must be logical scalar (true or false).');
end

% Check if input is valid
if strcmp(AQN.storage_mode,'limitedMemory') ~= 1
  error('AggQN: tryNewton incompatible with storage mode.  For aggregation, storage mode needs to be ''limitedMemory''.');
end

% Print message
if AQN.tryNewton ~= value && AQN.verbosity >= 1
  fprintf('AggQN: tryNewton set to %d\n',value);
end

% If all is OK, then set option
AQN.tryNewton = value;

end
