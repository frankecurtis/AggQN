% Copyright (C) 2021 Albert Berahas, Frank E. Curtis, Baoyu Zhou
%
% All Rights Reserved.
%
% Authors: Albert Berahas, Frank E. Curtis, Baoyu Zhou
%
% Method definition for AggQN class

function deleteDataLimitedMemory(AQN,index,reason)

% Delete pair
AQN.SY(:,index)  = [];
AQN.SY(index,:)  = [];
AQN.R(:,index)   = [];
AQN.R(index,:)   = [];
AQN.D(:,index)   = [];
AQN.D(index,:)   = [];
AQN.L(:,index)   = [];
AQN.L(index,:)   = [];
AQN.HS(:,index)  = [];
AQN.SHS(:,index) = [];
AQN.SHS(index,:) = [];
AQN.S(:,index)   = [];
AQN.Y(:,index)   = [];
AQN.rho(index)   = [];

% Print message
if AQN.verbosity >= 1
  fprintf('AggQN: Pair deleted, index %d; %s\n',index,reason);
end

end
