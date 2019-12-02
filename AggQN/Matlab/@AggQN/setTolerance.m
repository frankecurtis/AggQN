% Copyright (C) 2019 Albert Berahas, Frank E. Curtis, Baoyu Zhou
%
% All Rights Reserved.
%
% Authors: Albert Berahas, Frank E. Curtis, Baoyu Zhou
%
% Method definition for AggQN class

% Sets verbosity
function setTolerance(AQN,tol1,tol2,tol3,tol4,tol5,tol6)

% Check input type
if ~isnumeric([tol1,tol2,tol3,tol4,tol5]) || tol1 ~= abs(tol1) || tol2 ~= abs(tol2) || tol3 ~= abs(tol3) || tol4 ~= abs(tol4) || tol5 ~= abs(tol5) || tol6 ~= abs(tol6)
  error('AggQN: Invalid input to setTolerance.  Input must be nonnegative number.');
end

% If all is OK, then set option
AQN.cond_tol_1   = 10^tol1;
AQN.cond_tol_2   = 10^tol2;
AQN.cond_tol_3   = 10^tol3;
AQN.angle_tol_1  = tol4;
AQN.angle_tol_2  = tol5;
AQN.accuracy_tol = tol6;
AQN.diff_tol     = tol6;

end
