% Copyright (C) 2019 Albert Berahas, Frank E. Curtis, Baoyu Zhou
%
% All Rights Reserved.
%
% Authors: Albert Berahas, Frank E. Curtis, Baoyu Zhou
%
% Method definition for AggQN class

% Check whether accuracy reached
function flag = checkAccuracyNewton(AQN,vec)

if norm(vec,inf) > AQN.check_tol * max(1,AQN.startValue)
    flag = 0;
else
    flag = 1;
end

end