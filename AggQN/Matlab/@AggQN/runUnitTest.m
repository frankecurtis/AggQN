% Copyright (C) 2021 Albert Berahas, Frank E. Curtis, Baoyu Zhou
%
% All Rights Reserved.
%
% Authors: Albert Berahas, Frank E. Curtis, Baoyu Zhou
%
% Method definition for AggQN class

% Runs unit test
function runUnitTest(AQN,number,level,phi,phi_comb,N)

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
      ak = AQN.invM_j*[AQN.MA_j(1:k,k); zeros(size(AQN.S,2)-AQN.j-k+1,1)];

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
      fprintf('AggQN: Error in computation of a_{%4d,2}     : %e\n',size(AQN.Y,2)-AQN.j,abs(v'*AQN.invM_j*v - AQN.rhs_j(size(AQN.Y,2)-AQN.j,size(AQN.Y,2)-AQN.j)));
    end

    %%%%%%%%%%%%%%%%%%%%%%%%
    % Unit test for (3.30) %
    %%%%%%%%%%%%%%%%%%%%%%%%
  case 3

    % Print error
    if AQN.verbosity >= 1
      fprintf('AggQN: Error in computation of a_{%4d,2}^*   : %e\n',level,max(abs(phi(:,level:end)'*AQN.invM_j*phi_comb - AQN.rhs_j(level+1:end,level))));
    end

    %%%%%%%%%%%%%%%%%%%%%%%%
    % Unit test for (3.33) %
    %%%%%%%%%%%%%%%%%%%%%%%%
  case 4

    % Print error
    if AQN.verbosity >= 1
      fprintf('AggQN: Error in computation of a_{%4d,2}-bar : %e\n',level,max(abs(phi(:,level:end)'*AQN.invM_j*N)));
    end

  otherwise

    % Print bad input
    fprintf('AggQN: Oops... bad input to runUnitTest.\n');

end

end
