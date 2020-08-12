% Copyright (C) 2019 Albert Berahas, Frank E. Curtis, Baoyu Zhou
%
% All Rights Reserved.
%
% Authors: Albert Berahas, Frank E. Curtis, Baoyu Zhou
%
% Method definition for AggQN class

function addDataLimitedMemory(AQN,s,y,reason)

% Add to S*Y matrix
if size(AQN.S,2) >= 1
  AQN.SY = [AQN.SY (y'*AQN.S)'; s'*AQN.Y s'*y];
else
  AQN.SY = s'*y;
end

% Add to S and Y
AQN.S = [AQN.S s];
AQN.Y = [AQN.Y y];

% Add to rho
AQN.rho = [AQN.rho; 1/(AQN.SY(end,end))];

% Print message
if AQN.verbosity >= 1
  fprintf('AggQN: Pair added; %s\n',reason);
end

end