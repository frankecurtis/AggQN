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
msg = [msg newline '       AQN = AggQN(''limitedMemory'',initWv,initHv,n,m)'];
msg = [msg newline '             where storage mode ''limitedMemory'' means'];
msg = [msg newline '                   initWv is a handle of a function that computes matrix-vector products with the ''initial'' inverse Hessian,'];
msg = [msg newline '                   initHv is a handle of a function that computes matrix-vector products with the ''initial'' Hessian,'];
msg = [msg newline '                   n is the length of s and y vectors, and'];
msg = [msg newline '                   m is the (limited memory) history length'];
msg = [msg newline '       AQN = AggQN({''denseInverseHessian'',''denseHessian''},M)'];
msg = [msg newline '             where storage mode ''denseInverseHessian'' means'];
msg = [msg newline '                   M is initial inverse Hessian'];
msg = [msg newline '             where storage mode ''denseHessian'' means'];
msg = [msg newline '                   M is initial Hessian'];

end
