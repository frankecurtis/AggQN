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
msg = [msg newline];
msg = [msg newline '      AQN = AggQN(''SY'',I,m)'];
msg = [msg newline '            where storage mode ''SY'''];
msg = [msg newline '                  means initWv is a handle of a function that computes'];
msg = [msg newline '                  matrix-vector products with the "initial" inverse Hessian'];
msg = [msg newline '                  and initHv is a handle of a function that computes'];
msg = [msg newline '                  matrix-vector products with the "initial" Hessian'];
msg = [msg newline '                  and the history m must be a positive integer or inf'];
msg = [msg newline];
msg = [msg newline '      AQN = AggQN({''W'',''H''},M)'];
msg = [msg newline '            where storage mode ''W'''];
msg = [msg newline '                  means M is initial inverse Hessian'];
msg = [msg newline '            where storage mode ''H'''];
msg = [msg newline '                  means M is initial Hessian'];

end
