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
    
    % Compute Cholesky Factorization of matrix [AQN.S s]'*[AQN.S s]
    S_all = [AQN.S s]'*[AQN.S s];
    [AQN.L_SS_cond,~] = AQN.choleskyPerturb(S_all);
    eig_all = diag(AQN.L_SS_cond);
    cond_all = (max(eig_all)/min(eig_all))^2;
    
    % Add only if linearly independent and size of S still small
    if cond_all <= AQN.cond_tol_1
        
        AQN.j = 1;
        
        % Store pair data for aggregation
        AQN.s_j   = AQN.S(:,AQN.j);
        AQN.y_j   = AQN.Y(:,AQN.j);
        AQN.rho_j = AQN.rho(AQN.j);
        AQN.tau_j = S_all(2:end,2:end)\S_all(2:end,1);
        
        % Add data
        AQN.addDataSY(s,y);
          
        if (abs(AQN.s_j'*AQN.S(:,AQN.j+1:end)*AQN.tau_j)/(norm(AQN.s_j)*norm(AQN.S(:,AQN.j+1:end)*AQN.tau_j))) < cos(AQN.angle_tol_1) || size(AQN.S,2) <= 2
            
            if size(AQN.S,2) > AQN.m
                
                % Delete data for SY
                AQN.deleteDataSY(AQN.j);
                    
                % Set message
                msg = 'Sw1';
                
                % Set aggregation accuracy
                AQN.aggAccuracy = 0;
                
            else
                
                % Set message
                msg = 'Ad1';
                
                % Set aggregation accuracy
                AQN.aggAccuracy = 0;
                
                % Print message
                if AQN.verbosity >= 1
                    fprintf('AggQN: Pair added into inverse Hessian approximation.\n');
                end
                
            end
            
        else
            
            % Use _temp to store information
            AQN.S_temp = AQN.S;
            AQN.Y_temp = AQN.Y;
            AQN.SY_temp = AQN.SY;
            AQN.rho_temp = AQN.rho;
            AQN.L_SS_temp = AQN.L_SS;
            
            % Delete data
            AQN.deleteDataSY(AQN.j);
            
            % Check if we do a preconditioning
            if AQN.precondition == 1
                
                % Compute aggregation vector
                AQN.computeAggregationValues;
                
                % Get aggregation accuracy
                accuracy = AQN.getAggAccuracy;
                
                % Check do aggregation or not
                if accuracy < AQN.accuracy_tol
            
                    % Set new Y
                    AQN.Y(:,AQN.j:end-1) = AQN.Y(:,AQN.j:end-1) + AQN.HQ_j*AQN.A_j + AQN.y_j*AQN.b_j';
                    
                    % Set message
                    msg = 'Ag1';
                    
                else
                    
                    if size(AQN.S) >= AQN.m
                        
                        % Set message
                        msg = 'Sw2';
                        
                        % Set aggregation accuracy
                        AQN.aggAccuracy = 0;
                        
                    else
                        
                        % Set message
                        msg = 'Ad2';
                        
                        % Set aggregation accuracy
                        AQN.aggAccuracy = 0;
                    
                        % Get data back
                        AQN.S = AQN.S_temp;
                        AQN.Y = AQN.Y_temp;
                        AQN.SY = AQN.SY_temp;
                        AQN.rho = AQN.rho_temp;
                        AQN.L_SS = AQN.L_SS_temp;
                    end
                end

                
            else
            
                % Compute aggregation vector
                AQN.computeAggregationValuesOld;
                
                % Check do aggregation or not
                if AQN.checkFlag == 1
            
                    % Set new Y
                    AQN.Y(:,AQN.j:end-1) = AQN.Y(:,AQN.j:end-1) + AQN.HS_j*AQN.A_j + AQN.y_j*AQN.b_j';
                    
                    % Set message
                    msg = 'Ag1';
                    
                else
                    
                    % Compute aggregation vector
                    AQN.computeAggregationValuesWithNewton;
                    
                    % Check do aggregation or not
                    if AQN.checkFlag == 1
                    
                        % Set new Y
                        AQN.Y(:,AQN.j:end-1) = AQN.Y(:,AQN.j:end-1) + AQN.HS_j*AQN.A_j + AQN.y_j*AQN.b_j';
                    
                        % Set message
                        msg = 'AN1';
                        
                    else
                    
                        if size(AQN.S) >= AQN.m
                        
                            % Set message
                            msg = 'Sw2';
                        
                            % Set aggregation accuracy
                            AQN.aggAccuracy = 0;
                        
                        else
                        
                            % Set message
                            msg = 'Ad2';
                        
                            % Set aggregation accuracy
                            AQN.aggAccuracy = 0;
                    
                            % Get data back
                            AQN.S = AQN.S_temp;
                            AQN.Y = AQN.Y_temp;
                            AQN.SY = AQN.SY_temp;
                            AQN.rho = AQN.rho_temp;
                            AQN.L_SS = AQN.L_SS_temp;
                        end
                        
                    end
                end
            end
                       
            % Update S'*Y
            AQN.SY(:,AQN.j:end-1) = AQN.S'*AQN.Y(:,AQN.j:end-1);
            
            % Print message
            if AQN.verbosity >= 1
              fprintf('AggQN: Pair aggregated into inverse Hessian approximation.\n');
            end
            
        end   
      
    else
      
      %%%%%%%%%%%%%%
      % AGGREGATE! %
      %%%%%%%%%%%%%%

      % Throw away multiple (s,y) pairs if necessary
      recur = true;
      i1 = size(AQN.S,2);
            
      % To make more efficient to check less j.
      while recur
          
          % Set initial R_temp
          R_temp = sqrt(s'*s);
          
          for i = i1 : -1 : 1
              if i == 1
                  AQN.j = i;
                  i1 = i;
                  break;
              else        
                  R_temp = AQN.choleskyRowAdd(R_temp,[s AQN.S(:,end:-1:i)]'*AQN.S(:,i),size(R_temp,1)+1);
                  eig_temp = diag(R_temp);
                  cond_temp = (max(eig_temp)/min(eig_temp))^2;
                  if cond_temp > AQN.cond_tol_2
                      AQN.j = i;
                      i1 = i-1;
                      break;
                  end
              end
          end
          
          R_del = AQN.choleskyRowDel(AQN.L_SS_cond',AQN.j);
          eig_del = diag(R_del);
          cond_del = (max(eig_del)/min(eig_del))^2;
                
          if cond_del > AQN.cond_tol_3
              recur = true;
              AQN.deleteDataSY(AQN.j);
              AQN.L_SS_cond = (AQN.choleskyRowDel(AQN.L_SS_cond',AQN.j))';
          else
              recur = false;
          end
      end
            
      % Sanity check for being parallel with latest pair
      if AQN.j == size(AQN.S,2) && AQN.j == 1
                
          % Set message
          msg = 'Ad0';
          
          % Set aggregation accuracy
          AQN.aggAccuracy = 0;
                
          % Add data
          AQN.addDataSY(s,y);
                
      elseif AQN.j == size(AQN.S,2) 
                
          % Delete data
          AQN.deleteDataSY(size(AQN.S,2));
                
          % Add data
          AQN.addDataSY(s,y);
                
          % Set message
          msg = 'Pr2';
          
          % Set aggregation accuracy
          AQN.aggAccuracy = 0;
                
          % Print message
          if AQN.verbosity >= 1
              fprintf('AggQN: Pair aggregated into inverse Hessian approximation; parallel with previous pair based in cond_tol_2.\n');
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
          if (abs(AQN.s_j'*AQN.S(:,AQN.j+1:end)*AQN.tau_j)/(norm(AQN.s_j)*norm(AQN.S(:,AQN.j+1:end)*AQN.tau_j))) < cos(AQN.angle_tol_2)
            
            % Delete data for SY
            AQN.deleteDataSY(AQN.j);
                    
            % Set message
            msg = 'Swj';
            
            % Set aggregation accuracy
            AQN.aggAccuracy = 0;

            % Print message
            if AQN.verbosity >= 1
              fprintf('AggQN: Pair added into inverse Hessian approximation.\n');
            end
            
          else
            
            % Use _temp to store information
            AQN.S_temp = AQN.S;
            AQN.Y_temp = AQN.Y;
            AQN.SY_temp = AQN.SY;
            AQN.rho_temp = AQN.rho;
            AQN.L_SS_temp = AQN.L_SS; 
              
            % Delete data
            AQN.deleteDataSY(AQN.j);
            
            % Check if we do a preconditioning
            if AQN.precondition == 1
                
                % Compute aggregation vector
                AQN.computeAggregationValues;
                
                % Get aggregation accuracy
                accuracy = AQN.getAggAccuracy;
                
                 % Check do aggregation or not
                if accuracy < AQN.accuracy_tol
            
                    % Set new Y
                    AQN.Y(:,AQN.j:end-1) = AQN.Y(:,AQN.j:end-1) + AQN.HQ_j*AQN.A_j + AQN.y_j*AQN.b_j';
                    
                    % Set message
                    msg = 'Ag2';
                    
                else
                    
                    if size(AQN.S) >= AQN.m
                        
                        % Set message
                        msg = 'Sw3';
                        
                        % Set aggregation accuracy
                        AQN.aggAccuracy = 0;
                        
                    else
                        
                        % Set message
                        msg = 'Ad3';
                        
                        % Set aggregation accuracy
                        AQN.aggAccuracy = 0;
                    
                        % Get data back
                        AQN.S = AQN.S_temp;
                        AQN.Y = AQN.Y_temp;
                        AQN.SY = AQN.SY_temp;
                        AQN.rho = AQN.rho_temp;
                        AQN.L_SS = AQN.L_SS_temp;
                        
                    end
                end
                
            else
            
                % Compute aggregation vector
                AQN.computeAggregationValuesOld;

                % Check do aggregation or not
                if AQN.checkFlag == 1
            
                    % Set new Y
                    AQN.Y(:,AQN.j:end-1) = AQN.Y(:,AQN.j:end-1) + AQN.HS_j*AQN.A_j + AQN.y_j*AQN.b_j';
                    
                    % Set message
                    msg = 'Ag2';
                    
                else
                    
                    % Compute aggregation vector
                    AQN.computeAggregationValuesWithNewton;
                    
                    % Check do aggregation or not
                    if AQN.checkFlag == 1
                    
                        % Set new Y
                        AQN.Y(:,AQN.j:end-1) = AQN.Y(:,AQN.j:end-1) + AQN.HS_j*AQN.A_j + AQN.y_j*AQN.b_j';
                    
                        % Set message
                        msg = 'AN2';
                        
                    else
                    
                        if size(AQN.S) >= AQN.m
                        
                            % Set message
                            msg = 'Sw3';
                        
                            % Set aggregation accuracy
                            AQN.aggAccuracy = 0;
                        
                        else
                        
                            % Set message
                            msg = 'Ad3';
                        
                            % Set aggregation accuracy
                            AQN.aggAccuracy = 0;
                    
                            % Get data back
                            AQN.S = AQN.S_temp;
                            AQN.Y = AQN.Y_temp;
                            AQN.SY = AQN.SY_temp;
                            AQN.rho = AQN.rho_temp;
                            AQN.L_SS = AQN.L_SS_temp;
                        end
                        
                    end
                                
                end
                
            end
                       
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

