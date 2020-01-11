% Copyright (C) 2019 Albert Berahas, Frank E. Curtis, Baoyu Zhou
%
% All Rights Reserved.
%
% Authors: Albert Berahas, Frank E. Curtis, Baoyu Zhou
%
% Method definition for AggQN class

% Sets verbosity
function setPrecondition(AQN,value)

% Check input type
if length(value) ~= 1 || ~islogical(value)
  error('AggQN: Invalid input to setPrecondition.  Input must be logical scalar (true or false).');
end

% Check if input is valid
if strcmp(AQN.storage_mode,'SY') ~= 1
  error('AggQN: Aggregation incompatible with storage mode.  For aggregation, storage mode needs to be SY.');
end

% Print message
if AQN.precondition ~= value && AQN.verbosity >= 1
  fprintf('AggQN: Precondition set to %d\n',value);
end

% If all is OK, then set option
AQN.precondition = value;

end