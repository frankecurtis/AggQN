% Copyright (C) 2019 Albert Berahas, Frank E. Curtis, Baoyu Zhou
%
% All Rights Reserved.
%
% Authors: Albert Berahas, Frank E. Curtis, Baoyu Zhou
%
% Method definition for AggQN class

% Add pair, storage mode 'limitedMemory'
function msg = addPairLimitedMemory(AQN,s,y)

% Check aggregation option
if ~AQN.aggregate
  
  %%%%%%%%%%%%%%%%%%
  % NO AGGREGATION %
  %%%%%%%%%%%%%%%%%%
  
  % Add data
  AQN.addDataLimitedMemory(s,y,'aggregation is off');
  
  % Set message
  msg = 'Add';
  
  % Check size
  if size(AQN.S,2) > AQN.m
    
    % Delete data
    AQN.deleteDataLimitedMemory(1,'history limit reached');
    
    % Update message
    msg = 'Swp';
    
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
    AQN.addDataLimitedMemory(s,y,'first pair');
    
    % Set message
    msg = 'Add';
    
  else
    
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    % CHECK WHETHER TO AGGREGATE %
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    
    % Add data
    AQN.addDataLimitedMemory(s,y,'aggregation is on');
    
    % Set message
    msg = 'Add';
    
    % Loop backwards, checking whether to aggregate
    for index = size(AQN.S,2)-1:-1:1
      
      % Compute potential tau
      tau = (AQN.S(:,index+1:end)'*AQN.S(:,index+1:end))\(AQN.S(:,index+1:end)'*AQN.S(:,index));
      
      % Compute projection
      projection = AQN.S(:,index+1:end)*tau;
      
      % Compute error vector
      error = AQN.S(:,index) - projection;
      
      % Compute bool for aggregation
      aggregation = (norm(error) <= AQN.proj_tol * norm(projection));
      
      % Check whether to try aggregation
      if aggregation || (size(AQN.S,2) > AQN.m && index == 1 && AQN.Y(:,1)'*projection >= AQN.curv_tol*norm(AQN.Y(:,1))*norm(projection) && norm(error) <= AQN.proj_tol_loose * norm(projection))
        
        % Store data, potentially to restore it
        S_temp   = AQN.S;
        Y_temp   = AQN.Y;
        rho_temp = AQN.rho;
        SY_temp  = AQN.SY;
        
        % Set aggregation values
        AQN.j     = index;
        AQN.s_j   = AQN.S(:,AQN.j);
        AQN.y_j   = AQN.Y(:,AQN.j);
        AQN.rho_j = AQN.rho(AQN.j);
        AQN.tau_j = tau;
        
        % Delete data
        AQN.deleteDataLimitedMemory(AQN.j,'trying aggregation');
        
        % Check if parallel to last
        if AQN.j == size(AQN.S,2)
          
          % Set message
          msg = 'Par';
          
          % Increment aggregation counter
          AQN.count = AQN.count + 1;
          
        else
          
          % Initialize aggregation values
          AQN.computeAggregationValuesInitial;
          
          % Compute aggregation values "exactly"
          AQN.computeAggregationValuesExact;
          
          % Compute error
          [errorMaxAbs,~,~,~] = AQN.computeAggregationValuesError(AQN.A_j);
          
          % Check do aggregation or not
          if errorMaxAbs <= AQN.acc_tol
            
            % Save error
            AQN.error = errorMaxAbs;
            
            % Set new Y
            if ~AQN.onlySwap
              AQN.Y(:,AQN.j:end-1) = AQN.Y(:,AQN.j:end-1) + AQN.HS_j*AQN.A_j + AQN.y_j*AQN.b_j';
            end
            
            % Update S'*Y
            AQN.SY(:,AQN.j:end-1) = AQN.S'*AQN.Y(:,AQN.j:end-1);
            
            % Set message
            msg = 'Agg';
            
            % Increment aggrecation counter
            AQN.count = AQN.count + 1;
            
          else
            
            % Try Newton?
            if AQN.tryNewton
              
              % Compute aggregation values with Newton's method
              AQN.computeAggregationValuesNewton;
              
              % Compute error
              [errorMaxAbs,~,~,~] = AQN.computeAggregationValuesError(AQN.A_j);
              
            end
            
            % Check do aggregation or not
            if errorMaxAbs <= AQN.acc_tol || (AQN.tryNewton && errorMaxAbs <= AQN.acc_tol_newton)
              
              % Save error
              AQN.error = errorMaxAbs;
              
              % Set new Y
              if ~AQN.onlySwap
                AQN.Y(:,AQN.j:end-1) = AQN.Y(:,AQN.j:end-1) + AQN.HS_j*AQN.A_j + AQN.y_j*AQN.b_j';
              end
              
              % Update S'*Y
              AQN.SY(:,AQN.j:end-1) = AQN.S'*AQN.Y(:,AQN.j:end-1);
              
              % Set message
              msg = 'AgN';
              
              % Increment aggrecation counter
              AQN.count = AQN.count + 1;
              
            else
              
              % Aggregation failed, restore data
              AQN.S   = S_temp;
              AQN.Y   = Y_temp;
              AQN.rho = rho_temp;
              AQN.SY  = SY_temp;
              
            end
            
          end
          
        end
        
      end
      
      % Check message
      if strcmp(msg,'Agg') == 1
        break;
      end
      
    end
    
    % Check size
    if size(AQN.S,2) > AQN.m
      
      % Delete data
      AQN.deleteDataLimitedMemory(1,'history limit reached');
      
      % Update message
      msg = 'Swp';
      
    end
    
  end
  
end

end