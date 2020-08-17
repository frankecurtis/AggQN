% Copyright (C) 2019 Albert Berahas, Frank E. Curtis, Baoyu Zhou
%
% All Rights Reserved.
%
% Authors: Albert Berahas, Frank E. Curtis, Baoyu Zhou
%
% Method definition for AggQN class

function addDataLimitedMemory(AQN,s,y,reason)

% Add to compact form matrices
if size(AQN.S,2) >= 1
  AQN.SY  = [AQN.SY (y'*AQN.S)'; s'*AQN.Y s'*y];
  AQN.R   = triu(AQN.SY);
  AQN.D   = diag(diag(AQN.SY));
  AQN.L   = tril(AQN.SY,-1);
  AQN.HS  = [AQN.HS AQN.initHv(s)];
  SHs     = (s'*AQN.HS(:,1:end-1))';
  AQN.SHS = [AQN.SHS SHs; SHs' s'*AQN.HS(:,end)];
else
  AQN.SY  = s'*y;
  AQN.R   = AQN.SY;
  AQN.D   = AQN.SY;
  AQN.L   = 0;
  AQN.HS  = AQN.initHv(s);
  AQN.SHS = s'*AQN.HS;
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