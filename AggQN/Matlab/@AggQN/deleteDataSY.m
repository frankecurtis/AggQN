% Copyright (C) 2019 Albert Berahas, Frank E. Curtis, Baoyu Zhou
%
% All Rights Reserved.
%
% Authors: Albert Berahas, Frank E. Curtis, Baoyu Zhou
%
% Method definition for AggQN class

function deleteDataSY(AQN,index,reason)

% Delete pair
AQN.SY(:,index) = [];
AQN.SY(index,:) = [];
AQN.S(:,index)  = [];
AQN.Y(:,index)  = [];
AQN.rho(index)  = [];

% Print message
if AQN.verbosity >= 1
  fprintf('AggQN: Pair deleted, index %d; %s\n',index,reason);
end

end