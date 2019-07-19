% Copyright (C) 2019 Albert Berahas, Frank E. Curtis, Baoyu Zhou
%
% All Rights Reserved.
%
% Authors: Albert Berahas, Frank E. Curtis, Baoyu Zhou
%
% Method definition for AggQN class

% Add pair, storage mode 'SY'
function msg = addPairSY(AQN,s,y)

% Check aggregation option
if ~AQN.aggregate
  
  %%%%%%%%%%%%%%%%%%
  % NO AGGREGATION %
  %%%%%%%%%%%%%%%%%%
  
  % Set message
  msg = 'Add';
  
  % Check size
  if size(AQN.S,2) >= AQN.m
    
    % Delete data
    AQN.deleteDataSY(1);
    
    % Update message
    msg = 'Swp';
    
    % Print message
    if AQN.verbosity >= 1
      fprintf('AggQN: Pair deleted from inverse Hessian approximation; history limit reached.\n');
    end
    
  end
    
  % Add data
  AQN.addDataSY(s,y);
    
  % Print message
  if AQN.verbosity >= 1
    fprintf('AggQN: Pair added into inverse Hessian approximation.\n');
  end
  
else
  
  %%%%%%%%%%%%%%%
  % AGGREGATION %
  %%%%%%%%%%%%%%%
  
  % Check cases
  if size(AQN.S,2) == 0
    
    %%%%%%%%%%%%%%%%
    % NO PAIRS YET %
    %%%%%%%%%%%%%%%%
    
    % Add data
    AQN.addDataSY(s,y);
    
    % Set message
    msg = 'Add';
    
    % Print message
    if AQN.verbosity >= 1
      fprintf('AggQN: Pair added into inverse Hessian approximation; first pair.\n');
    end
    
  elseif abs(AQN.S(:,end)'*s)/(norm(AQN.S(:,end))*norm(s)) >= 1 - AQN.parallel_tol
    
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    % PARALLEL WITH MOST RECENT PAIR %
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    
    % Set actual pair aggregated
    AQN.j = size(AQN.S,2);
    
    % Delete data
    AQN.deleteDataSY(size(AQN.S,2));
    
    % Add data
    AQN.addDataSY(s,y);
    
    % Set message
    msg = 'Swp';
    
    % Print message
    if AQN.verbosity >= 1
      fprintf('AggQN: Pair aggregated into inverse Hessian approximation; parallel with previous pair based on parallel_tol.\n');
    end
    
  else
    
    %%%%%%%%%%%%%%%%%%%%
    % ADD OR AGGREGATE %
    %%%%%%%%%%%%%%%%%%%%
    
    % In the following two steps,
    % values are "flipped", then "flipped" back,
    % due to reverse order in factorization of S'*S
    
    % Compute rhs = S'*s
    %rhs = (s'*AQN.S)';
    %rhs = flipud(rhs);
    
    % Solve S'*S*v = L*L'*v = S'*s
    %  i.e.          L*   u = S'*s
    %  then            L'*v = u
    %u = AQN.L_SS\rhs;
    %v = ((u')/AQN.L_SS)';
    %v = flipud(v);
    
    % Add only if linearly independent and size of S still small
    %if norm(AQN.S*v - s)/norm(AQN.S*v) > AQN.lin_ind_tol && size(AQN.S,2) < min(AQN.m,AQN.n)
    if cond([AQN.S s]'*[AQN.S s]) <= AQN.cond_tol_1 && size(AQN.S,2) < AQN.m
      
      %%%%%%%%%%%%
      % ADD ONLY %
      %%%%%%%%%%%%
      
      % Add data
      AQN.addDataSY(s,y);
      
      % Set message
      msg = 'Add';
      
      % Print message
      if AQN.verbosity >= 1
        fprintf('AggQN: Pair added into inverse Hessian approximation.\n');
      end
      
    else
      
      %%%%%%%%%%%%%%
      % AGGREGATE! %
      %%%%%%%%%%%%%%
      
      % Set bool for doing rotation
      % ... true if steps are linearly independent, false otherwise
      %do_rotation = (norm(AQN.S*v - s)/norm(AQN.S*v) > AQN.lin_ind_tol);
      
      % Find pair to aggregate (TO DO: MAKE MORE EFFICIENT)
      for i = 1:size(AQN.S,2)
        if i == size(AQN.S,2) || cond([AQN.S(:,end-i+1:end) s]'*[AQN.S(:,end-i+1:end) s]) > AQN.cond_tol_2
          AQN.j = size(AQN.S,2) - i + 1;
          break;
        end
      end
      
      % Sanity check for being parallel with latest pair
      if AQN.j == size(AQN.S,2)
        
        % Delete data
        AQN.deleteDataSY(size(AQN.S,2));
        
        % Add data
        AQN.addDataSY(s,y);
        
        % Set message
        msg = 'Swp';
        
        % Print message
        if AQN.verbosity >= 1
          fprintf('AggQN: Pair aggregated into inverse Hessian approximation; parallel with previous pair based in cond_tol_2.\n');
        end
        
      else
        
        % Check condition number of potential S
        if cond([AQN.S(:,1:AQN.j-1) AQN.S(:,AQN.j+1:end)]'*[AQN.S(:,1:AQN.j-1) AQN.S(:,AQN.j+1:end)]) > AQN.cond_tol_3
          
          % Set message
          msg = 'Adj';
          
          % Check history length
          if size(AQN.S,2) >= AQN.m
            
            % Delete oldest pair
            AQN.deleteDataSY(1);
            
            % Update message
            msg = 'Swj';
          
          end
          
          % Add new data
          AQN.addDataSY(s,y);
          
          % Print message
          if AQN.verbosity >= 1
            fprintf('AggQN: Pair swapped into inverse Hessian approximation; oldest pair forgotten.\n');
          end
          
        else
                    
          % Store pair data for aggregation
          AQN.s_j   = AQN.S(:,AQN.j);
          AQN.y_j   = AQN.Y(:,AQN.j);
          AQN.rho_j = AQN.rho(AQN.j);
          
          % Add data
          AQN.addDataSY(s,y);
          
          % In the following two steps,
          % values are "flipped", then "flipped" back,
          % due to reverse order in factorization of S'*S
          
          % Compute rhs = S'*s
          rhs = (AQN.s_j'*AQN.S(:,AQN.j+1:end))';
          rhs = flipud(rhs);
          
          % Solve S'*S*tau = L*L'*tau  = S'*s
          %  i.e.            L*   temp = S'*s
          %  then              L'*tau  = temp
          temp      = AQN.L_SS(1:size(AQN.S,2)-AQN.j,1:size(AQN.S,2)-AQN.j)\rhs;
          AQN.tau_j = ((temp')/AQN.L_SS(1:size(AQN.S,2)-AQN.j,1:size(AQN.S,2)-AQN.j))';
          AQN.tau_j = flipud(AQN.tau_j);
          
          % Check angle
          if 0 && (abs(AQN.s_j'*AQN.S(:,AQN.j+1:end)*AQN.tau_j)/(norm(AQN.s_j)*norm(AQN.S(:,AQN.j+1:end)*AQN.tau_j))) < 0.75
            
            % Increase history length
            AQN.m = AQN.m + 1;
            
            % Set message
            msg = 'Ad2';

            % Print message
            if AQN.verbosity >= 1
              fprintf('AggQN: Pair added into inverse Hessian approximation.\n');
            end
            
          else
            
            % Delete data
            AQN.deleteDataSY(AQN.j);
                        
            % Set message
            msg = 'Agg';
            
%             % Set projected then rotated s_j
%             s_j_rotated = AQN.S(:,AQN.j:end)*AQN.tau_j;
%             AQN.tau_j = AQN.tau_j * norm(AQN.s_j) / norm(s_j_rotated);
%             s_j_rotated = s_j_rotated * norm(AQN.s_j) / norm(s_j_rotated);          
%             
%             % Set rotation matrix
%             outer_product   = AQN.s_j*s_j_rotated';
%             [V,D]           = eig(outer_product' + outer_product);
%             [~,index]       = min(diag(D));
%             eigen_vector    = V(:,index);
%             rotation_matrix = eye(length(AQN.s_j)) - 2*(eigen_vector*eigen_vector');
%             
% %             % Only for testing procondition testcase!
% %             rotation_matrix = eye(length(AQN.s_j));
%             
%             
%             % Rotate vectors
%             AQN.s_j = s_j_rotated;
%             AQN.y_j = rotation_matrix * AQN.y_j;
            
            % Set initial Hessian inverse
            W1 = AQN.initWv(eye(AQN.n));
            W2 = W1;
            
            % Set true rotation S,Y pairs
            if AQN.j == 1
                S1 = [AQN.s_j AQN.S];
                Y1 = [AQN.y_j AQN.Y];
            else
                S1 = [AQN.S(:,1:AQN.j-1) AQN.s_j AQN.S(:,AQN.j:end)];
                Y1 = [AQN.Y(:,1:AQN.j-1) AQN.y_j AQN.Y(:,AQN.j:end)];
            end
            
            if AQN.precondition == 1
                % Compute aggregation vector
                AQN.computeAggregationValues;
            
                % Set new Y
                AQN.Y(:,AQN.j:end-1) = AQN.Y(:,AQN.j:end-1) + AQN.HQ_j*AQN.A_j + AQN.y_j*AQN.b_j';
                
            else
            
                % Compute aggregation vector
                AQN.computeAggregationValuesOld;
            
                % Set new Y
                AQN.Y(:,AQN.j:end-1) = AQN.Y(:,AQN.j:end-1) + AQN.HS_j*AQN.A_j + AQN.y_j*AQN.b_j';
            end
             
            % Set S,Y pairs after aggregation
            S2 = AQN.S;
            Y2 = AQN.Y;
            
            % compare Hessian updates by (S1,Y1) and (S2,Y2)
            for i = 1:size(S1,2)
                U1 = eye(AQN.n) - (Y1(:,i)*S1(:,i)')/(Y1(:,i)'*S1(:,i));
                V1 = (S1(:,i)*S1(:,i)')/(Y1(:,i)'*S1(:,i));
                W1 = U1'*W1*U1 + V1;
            end
            
            for i = 1:size(S2,2)
                U2 = eye(AQN.n) - (Y2(:,i)*S2(:,i)')/(Y2(:,i)'*S2(:,i));
                V2 = (S2(:,i)*S2(:,i)')/(Y2(:,i)'*S2(:,i));
                W2 = U2'*W2*U2 + V2;
            end
            
            AQN.sum_diff = [AQN.sum_diff min(norm(W1-W2,inf)/norm(W1,inf),norm(W1-W2,inf))];
            
            % Update S'*Y
            AQN.SY(:,AQN.j:end-1) = AQN.S'*AQN.Y(:,AQN.j:end-1);
            
            % Print message
            if AQN.verbosity >= 1
              fprintf('AggQN: Pair aggregated into inverse Hessian approximation.\n');
            end
            
          end
          
        end
        
      end
      
    end
    
  end
  
end

end
