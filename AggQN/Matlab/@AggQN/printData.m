% Copyright (C) 2021 Albert Berahas, Frank E. Curtis, Baoyu Zhou
%
% All Rights Reserved.
%
% Authors: Albert Berahas, Frank E. Curtis, Baoyu Zhou
%
% Method definition for AggQN class

% Print data (for debugging)
function printData(AQN)

% Print scalar data
fprintf('AggQN: Printing data...\n');
fprintf('Storage mode = %s\n',AQN.storage_mode);
fprintf('Aggregation  = %s\n',mat2str(AQN.aggregate));
fprintf('Debug        = %s\n',mat2str(AQN.debug));
fprintf('OnlySwap     = %s\n',mat2str(AQN.onlySwap));
fprintf('TryNewton    = %s\n',mat2str(AQN.tryNewton));
fprintf('Verbosity    = %d\n',AQN.verbosity);
fprintf('History      = %d\n',AQN.m);
fprintf('Size         = %d\n',AQN.n);

% Print (inverse) Hessian data
if strcmp(AQN.storage_mode,'limitedMemory') == 1
  if size(AQN.S,2) >= 1
    fprintf('S            = \n');
    disp(AQN.S);
    fprintf('Y            = \n');
    disp(AQN.Y);
    fprintf('rho          = \n');
    disp(AQN.rho);
    fprintf('SY           = \n');
    disp(AQN.SY);
    rho = zeros(size(AQN.S,2),1);
    for i = 1:size(AQN.S,2)
      rho(i) = 1/(AQN.S(:,i)'*AQN.Y(:,i));
    end
    fprintf('Error in rho                  = %e\n',max(abs(AQN.rho - rho)));
    SY = zeros(size(AQN.S,2),size(AQN.Y,2));
    for i = 1:size(AQN.S,2)
      for j = 1:size(AQN.Y,2)
        SY(i,j) = AQN.S(:,i)'*AQN.Y(:,j);
      end
    end
    fprintf('Error in S''*  Y               = %e\n',max(max(abs(AQN.SY - SY))));
  end
end
if strcmp(AQN.storage_mode,'denseInverseHessian') == 1
  fprintf('W            = \n');
  disp(AQN.W);
end
if strcmp(AQN.storage_mode,'denseHessian') == 1
  fprintf('H            = \n');
  disp(AQN.H);
end

end
