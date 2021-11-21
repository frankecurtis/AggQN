% Copyright (C) 2021 Albert Berahas, Frank E. Curtis, Baoyu Zhou
%
% All Rights Reserved.
%
% Authors: Albert Berahas, Frank E. Curtis, Baoyu Zhou
%
% Method definition for AggQN class

% Sets aggregation
function setOnlySwap(AQN,value)

% Check input type
if length(value) ~= 1 || ~islogical(value)
  error('AggQN: Invalid input to setOnlyRemove.  Input must be logical scalar (true or false).');
end

% Check if input is valid
if strcmp(AQN.storage_mode,'limitedMemory') ~= 1
  error('AggQN: onlySwap incompatible with storage mode.  For onlySwap, storage mode needs to be ''limitedMemory''.');
end

% Print message
if AQN.onlySwap ~= value && AQN.verbosity >= 1
  fprintf('AggQN: onlySwap set to %d\n',value);
end

% If all is OK, then set option
AQN.onlySwap = value;

end
