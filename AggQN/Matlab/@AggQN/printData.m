% Copyright (C) 2019 Albert Berahas, Frank E. Curtis, Baoyu Zhou
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
fprintf('Adaptivity   = %s\n',mat2str(AQN.adaptive));
fprintf('Aggregation  = %s\n',mat2str(AQN.aggregate));
fprintf('Verbosity    = %d\n',AQN.verbosity);
fprintf('History      = %d\n',AQN.m);
fprintf('Size         = %d\n',AQN.n);

% Print (inverse) Hessian data
if strcmp(AQN.storage_mode,'SY') == 1
  if size(AQN.S,2) >= 1
    fprintf('S            = \n');
    disp(AQN.S);
    fprintf('Y            = \n');
    disp(AQN.Y);
    fprintf('rho          = \n');
    disp(AQN.rho);
    fprintf('SY           = \n');
    disp(AQN.SY);
    fprintf('L_SS         = \n');
    disp(AQN.L_SS);
    fprintf('L_SHS        = \n');
    disp(AQN.L_SHS);
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
    S_reverse = zeros(size(AQN.S));
    for i = 1:size(AQN.S,2)
      S_reverse(:,i) = AQN.S(:,end+1-i);
    end
    fprintf('Error in S''*  S factorization = %e\n',max(max(abs(S_reverse'*S_reverse     - AQN.L_SS *AQN.L_SS' ))));
    fprintf('Error in S''*H*S factorization = %e\n',max(max(abs(AQN.S'*AQN.initHv(AQN.S) - AQN.L_SHS*AQN.L_SHS'))));
  end
end
if strcmp(AQN.storage_mode,'W') == 1
  fprintf('W            = \n');
  disp(AQN.W);
end
if strcmp(AQN.storage_mode,'H') == 1
  fprintf('H            = \n');
  disp(AQN.H);
end

end
