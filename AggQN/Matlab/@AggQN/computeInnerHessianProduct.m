% Copyright (C) 2021 Albert Berahas, Frank E. Curtis, Baoyu Zhou
%
% All Rights Reserved.
%
% Authors: Albert Berahas, Frank E. Curtis, Baoyu Zhou
%
% Method definition for AggQN class

% Compute product with 'H_j'
function Hjv = computeInnerHessianProduct(AQN,v)

% Compute compact form product
temp = [AQN.SHS(1:AQN.j-1,1:AQN.j-1) AQN.L(1:AQN.j-1,1:AQN.j-1); AQN.L(1:AQN.j-1,1:AQN.j-1)' -AQN.D(1:AQN.j-1,1:AQN.j-1)]\([AQN.HS(:,1:AQN.j-1)'; AQN.Y(:,1:AQN.j-1)']*v);

% Compute product
Hjv = AQN.initHv(v) - [AQN.HS(:,1:AQN.j-1) AQN.Y(:,1:AQN.j-1)]*temp;

end
