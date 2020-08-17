% Copyright (C) 2019 Albert Berahas, Frank E. Curtis, Baoyu Zhou
%
% All Rights Reserved.
%
% Authors: Albert Berahas, Frank E. Curtis, Baoyu Zhou
%
% Method definition for AggQN class

% Compute product with 'W_j'
function Wjv = computeInnerInverseHessianProduct(AQN,v)

% Construct initial products
Wv  = AQN.initWv(v);
Stv = (v'*AQN.S(:,1:AQN.j-1))';

% Construct intermediate products
RinvStv           = AQN.R(1:AQN.j-1,1:AQN.j-1)\Stv;
YRinvStv          = AQN.Y(:,1:AQN.j-1)*RinvStv;
WYRinvStv         = AQN.initWv(YRinvStv);
YtWv              = (Wv'*AQN.Y(:,1:AQN.j-1))';
RinvtYtWv         = (YtWv'/AQN.R(1:AQN.j-1,1:AQN.j-1))';
SRinvtYtWv        = AQN.S(:,1:AQN.j-1)*RinvtYtWv;
DRinvStv          = AQN.D(1:AQN.j-1,1:AQN.j-1)*RinvStv;
RinvtDRinvStv     = (DRinvStv'/AQN.R(1:AQN.j-1,1:AQN.j-1))';
SRinvtDRinvStv    = AQN.S(:,1:AQN.j-1)*RinvtDRinvStv;
YtWYRinvStv       = (WYRinvStv'*AQN.Y(:,1:AQN.j-1))';
RinvtYtWYRinvStv  = (YtWYRinvStv'/AQN.R(1:AQN.j-1,1:AQN.j-1))';
SRinvtYtWYRinvStv = AQN.S(:,1:AQN.j-1)*RinvtYtWYRinvStv;

% Compute product
Wjv = Wv - WYRinvStv - SRinvtYtWv + SRinvtDRinvStv + SRinvtYtWYRinvStv;

end