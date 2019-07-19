% Copyright (C) 2019 Albert Berahas, Frank E. Curtis, Baoyu Zhou
%
% All Rights Reserved.
%
% Authors: Albert Berahas, Frank E. Curtis, Baoyu Zhou
%
% Method definition for AggQN class

% Inverse Hessian
function W = inverseHessian(AQN)

% Check storage mode
if strcmp(AQN.storage_mode,'SY') == 1
  
  % TO DO: Make more efficient in terms of storage and computational costs
  
  % Check size
  if size(AQN.S,2) >= 1
    
    % Construct R and D
    R = zeros(size(AQN.S,2),size(AQN.S,2));
    for i = 1:size(R,1)
      for j = i:size(R,2)
        R(i,j) = AQN.S(:,i)'*AQN.Y(:,j);
      end
    end
    D = diag(diag(R));
    
    % Construct intermediate matrices
    RinvSt   = R\(AQN.S');
    YRinvSt  = AQN.Y*RinvSt;
    WYRinvSt = AQN.initWv(eye(AQN.n))*YRinvSt;
    
    % Construct inverse Hessian
    W = AQN.initWv(eye(AQN.n)) - WYRinvSt - WYRinvSt' + RinvSt'*D*RinvSt + YRinvSt'*WYRinvSt;
    
    % Symmetrize
    W = (W + W')/2;
    
  else
    
    % Construct inverse Hessian
    W = AQN.initWv(eye(AQN.n));
    
  end
  
  % Print message
  warning('AggQN: Constructing inverse Hessian for storage mode ''SY'' is inefficient!  Consider computing inverse-Hessian-products instead.');
  
elseif strcmp(AQN.storage_mode,'W') == 1
  
  % Set return value
  W = (AQN.W + AQN.W')/2;
  
else % strcmo(AQN.storage_mode,'H') == 1
  
  % Set return value
  W = inv((AQN.H + AQN.H')/2);
  
  % Print message
  warning('AggQN: Constructing inverse Hessian for storage mode ''H'' is inefficient and inaccurate!');
  
end

end