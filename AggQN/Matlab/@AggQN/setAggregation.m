% Copyright (C) 2021 Albert Berahas, Frank E. Curtis, Baoyu Zhou
%
% All Rights Reserved.
%
% Authors: Albert Berahas, Frank E. Curtis, Baoyu Zhou
%
% Method definition for AggQN class

% Sets aggregation
function setAggregation(AQN,value)

% Check input type
if length(value) ~= 1 || ~islogical(value)
  error('AggQN: Invalid input to setAggregation.  Input must be logical scalar (true or false).');
end

% Check if input is valid
if strcmp(AQN.storage_mode,'limitedMemory') ~= 1
  error('AggQN: Aggregation incompatible with storage mode.  For aggregation, storage mode needs to be ''limitedMemory''.');
end

% Check if pairs have been initialized
if strcmp(AQN.storage_mode,'limitedMemory') == 1 && ~isempty(AQN.S)
  error('AggQN: Aggregation cannot be changed after a pair has already been added.');
end

% Print message
if AQN.aggregate ~= value && AQN.verbosity >= 1
  fprintf('AggQN: Aggregation set to %d\n',value);
end

% If all is OK, then set option
AQN.aggregate = value;

end
