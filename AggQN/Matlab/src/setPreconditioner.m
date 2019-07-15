% Copyright (C) 2019 Albert Berahas, Frank E. Curtis, Baoyu Zhou
%
% All Rights Reserved.
%
% Authors: Albert Berahas, Frank E. Curtis, Baoyu Zhou
%
% Method definition for AggQN class

% Sets verbosity
function setPreconditioner(AQN,value)

% Check input type
if length(value) ~= 1 || ~isnumeric(value) || value ~= abs(value)
  error('AggQN: Invalid input to setVerbosity.  Input must be nonnegative integer.');
end

% If all is OK, then set option
AQN.precondition = value;

end