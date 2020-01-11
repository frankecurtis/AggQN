% Copyright (C) 2019 Albert Berahas, Frank E. Curtis, Baoyu Zhou
%
% All Rights Reserved.
%
% Authors: Albert Berahas, Frank E. Curtis, Baoyu Zhou
%
% Method definition for AggQN class

% Sets adaptivity
function setAdaptivity(AQN,value)

% Check input type
if length(value) ~= 1 || ~islogical(value)
  error('AggQN: Invalid input to setAdaptivity.  Input must be logical scalar (true or false).');
end

% Check if input is valid
if strcmp(AQN.storage_mode,'SY') ~= 1
  error('AggQN: Adaptivity incompatible with storage mode.  For adaptivity, storage mode needs to be SY.');
end

% Check if pairs have been initialized
if strcmp(AQN.storage_mode,'SY') == 1 && ~isempty(AQN.S)
  error('AggQN: Adaptivity cannot be changed after a pair has already been added.');
end

% Print message
if AQN.adaptive ~= value && AQN.verbosity >= 1
  fprintf('AggQN: Adaptivity set to %d\n',value);
end

% If all is OK, then set option
AQN.adaptive = value;

end