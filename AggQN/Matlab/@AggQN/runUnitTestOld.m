% Copyright (C) 2019 Albert Berahas, Frank E. Curtis, Baoyu Zhou
%
% All Rights Reserved.
%
% Authors: Albert Berahas, Frank E. Curtis, Baoyu Zhou
%
% Method definition for AggQN class

% Runs unit test
function runUnitTestOld(AQN,number,invM,vec_temp,rhs,level,phi,phi_comb,N)

% Check whether in debug mode
if AQN.debug == false, return; end

% Switch on number
switch number

  %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
  % Unit test for "top part" of Q*A, which %
  % is supposed to satisfy (3.13)          %
  %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
  case 1
    
    % Loop over indices
    for k = 1:size(AQN.Y,2)-AQN.j
      
      % Set temporary vector
      ak = invM*[AQN.MA_j(1:k,k); zeros(size(AQN.S,2)-AQN.j-k+1,1)];
      
      % Print error
      if AQN.verbosity >= 1
        fprintf('AggQN: Error in computation of a_{%4d,1}     : %e\n',k,norm(AQN.S(:,AQN.j:AQN.j+k-1)'*AQN.HS_j*ak + AQN.b_j(k)*AQN.S(:,AQN.j:AQN.j+k-1)'*AQN.y_j));
      end
      
    end
    
  %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
  % Unit test for "first quadratic"     %
  % corresponding to l = m-1 in (3.21c) %
  %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
  case 2
    
    % Set temporary vector
    v = [zeros(size(AQN.Y,2)-AQN.j,1);
         AQN.MA_j(size(AQN.Y,2)-AQN.j+1,size(AQN.Y,2)-AQN.j) + AQN.S(:,end)'*(AQN.b_j(size(AQN.Y,2)-AQN.j)*AQN.y_j + AQN.Y(:,size(AQN.Y,2)-1))];
    
    % Print error
    if AQN.verbosity >= 1
      fprintf('AggQN: Error in computation of a_{%4d,2}     : %e\n',size(AQN.Y,2)-AQN.j,abs(v'*invM*v - rhs(size(AQN.Y,2)-AQN.j,size(AQN.Y,2)-AQN.j)));
    end
  
  %%%%%%%%%%%%%%%%%%%%%%%%
  % Unit test for (3.30) %
  %%%%%%%%%%%%%%%%%%%%%%%%
  case 3
    
    % Print error
    if AQN.verbosity >= 1
      fprintf('AggQN: Error in computation of a_{%4d,2}^*   : %e\n',level,max(abs(phi(:,level:end)'*invM*phi_comb - rhs(level+1:end,level))));
    end
  
  %%%%%%%%%%%%%%%%%%%%%%%%
  % Unit test for (3.33) %
  %%%%%%%%%%%%%%%%%%%%%%%%
  case 4
    
    % Print error
    if AQN.verbosity >= 1
      fprintf('AggQN: Error in computation of a_{%4d,2}-bar : %e\n',level,max(abs(phi(:,level:end)'*invM*N)));
    end
  
  %%%%%%%%%%%%%%%%%%%%%%%
  % Unit test for (3.9) %
  %%%%%%%%%%%%%%%%%%%%%%%
  case 5
    
%     % Construct H_{1:j-1} as AQN_H.hessian
%     H_matr = AQN.initHv(eye(AQN.n));
%     AQN_H  = AggQN('H',H_matr);
%     AQN_H.setVerbosity(0);
%     for i = 1:AQN.j-1
%       AQN_H.addPair(AQN.S(:,i),AQN.Y(:,i));
%     end
%     
%     % Compute Ytilde
%     Ytilde = AQN.Y;
%     Ytilde(:,AQN.j:end-1) = AQN.Y(:,AQN.j:end-1) + AQN_H.hessian*AQN.S(:,AQN.j:end)*AQN.A_j + AQN.y_j*AQN.b_j';
%     
%     % Compute "full" b
%     b = [zeros(AQN.j-1,1); AQN.b_j; 0];
%     
%     % Compute "full" A
%     A = [AQN.A_j zeros(size(AQN.S(:,AQN.j:end),2),1)];
%     
%     % Construct matrices
%     matrix1 = triu(AQN.SY(AQN.j:end,AQN.j:end)) - triu(AQN.S(:,AQN.j:end)'*Ytilde(:,AQN.j:end));
%     matrix2 = b(AQN.j:end) + AQN.rho_j*tril(AQN.SY(AQN.j:end,AQN.j:end),-1)'*AQN.tau_j;
%     matrix3 = abs((Ytilde(:,AQN.j:end) - AQN.Y(:,AQN.j:end))'*inv(AQN_H.hessian)*(Ytilde(:,AQN.j:end)-AQN.Y(:,AQN.j:end)) - ...
%               (1+AQN.rho_j*AQN.y_j'*inv(AQN_H.hessian)*AQN.y_j)/AQN.rho_j * (b(AQN.j:end)*b(AQN.j:end)') + ...
%               A'*tril(AQN.SY(AQN.j:end,AQN.j:end),-1) + ...
%               tril(AQN.SY(AQN.j:end,AQN.j:end),-1)'*A);
%     
%     % Print errors
%     if AQN.verbosity >= 1
%       fprintf('AggQN: Error in satisfaction of R = Rtilde    : %e\n',max(max(matrix1)));
%       fprintf('AggQN: Error in satisfaction of b definition  : %e\n',max(abs(matrix2)));
%       fprintf('AggQN: Error in satisfaction of quadratic     : %e\n',max(max(matrix3)));
%     end
    
    % Print errors
    if AQN.verbosity >= 1
      fprintf('AggQN: Error in satisfaction of quadratic-linear system    : %e\n',norm(vec_temp,inf));
    end

    AQN.aggAccuracy = norm(vec_temp,inf);
  
  otherwise
    
    % Print bad input
    fprintf('AggQN: Oops... bad input to runUnitTest.\n');
  
end

end