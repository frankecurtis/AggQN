% Copyright (C) 2019 Albert Berahas, Frank E. Curtis, Baoyu Zhou
%
% All Rights Reserved.
%
% Authors: Albert Berahas, Frank E. Curtis, Baoyu Zhou
%
% Method definition for AggQN class

function msg = usageMessage

% Set usage message
msg = 'AggQN: Usage: Constructors are';
msg = [msg];
msg = strcat(msg,'; ','AQN = AggQN(''SY'',I,m)');
msg = strcat(msg,'; ','where storage mode ''SY''');
msg = strcat(msg,'; ','means initWv is a handle of a function that computes');
msg = strcat(msg,'; ','matrix-vector products with the "initial" inverse Hessian');
msg = strcat(msg,'; ','and initHv is a handle of a function that computes');
msg = strcat(msg,'; ','matrix-vector products with the "initial" Hessian');
msg = strcat(msg,'; ','and the history m must be a positive integer or inf');
msg = [msg];
msg = strcat(msg,'; ','AQN = AggQN({''W'',''H''},M)');
msg = strcat(msg,'; ','where storage mode ''W''');
msg = strcat(msg,'; ','means M is initial inverse Hessian');
msg = strcat(msg,'; ','where storage mode ''H''');
msg = strcat(msg,'; ','means M is initial Hessian');

end
