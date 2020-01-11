% Copyright (C) 2019 Albert Berahas, Frank E. Curtis, Baoyu Zhou
%
% All Rights Reserved.
%
% Authors: Albert Berahas, Frank E. Curtis, Baoyu Zhou
%
% Method definition for AggQN class

function msg = usageMessage

% Set usage message
msg = 'AggQN: Usage:';
msg = [msg newline '       AQN = AggQN(''SY'',I,m)'];
msg = [msg newline '             where storage mode ''SY'' means'];
msg = [msg newline '                   initWv is a handle of a function that computes matrix-vector products with the ''initial'' inverse Hessian,'];
msg = [msg newline '                   initHv is a handle of a function that computes matrix-vector products with the ''initial'' Hessian,'];
msg = [msg newline '                   n is the length of s and y vectors, and'];
msg = [msg newline '                   m is the history length (for limited memory methods)'];
msg = [msg newline '       AQN = AggQN({''W'',''H''},M)'];
msg = [msg newline '             where storage mode ''W'' means'];
msg = [msg newline '                   M is initial inverse Hessian'];
msg = [msg newline '             where storage mode ''H'' means'];
msg = [msg newline '                   M is initial Hessian'];

end
