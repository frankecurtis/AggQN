% Copyright (C) 2019 Albert Berahas, Frank E. Curtis, Baoyu Zhou
%
% All Rights Reserved.
%
% Authors: Albert Berahas, Frank E. Curtis, Baoyu Zhou
%
% Method definition for AggQN class

% Compute aggregation values 'A' and 'b'
function computeAggregationValuesWithNewton(AQN)

% Construct H_{1:j-1} as AQN_H.hessian
H_matr = AQN.initHv(eye(AQN.n));
AQN_H  = AggQN('H',H_matr);
AQN_H.setVerbosity(0);
for i = 1:AQN.j-1
  AQN_H.addPair(AQN.S(:,i),AQN.Y(:,i));
end

% Construct HS_j
AQN.HS_j = AQN_H.hessian*AQN.S(:,AQN.j:end);

% Construct SHS_j
SHS_j = AQN.S(:,AQN.j:end)'*AQN.HS_j;

%%%%%%%%%%%%%
% Compute b %
%%%%%%%%%%%%%

% Compute b_j = -rho_j * (S(:,j:end)'*Y(:,j:end-1) - R(:,j:end-1))' * tau_j
AQN.b_j = -AQN.rho_j * tril(AQN.SY(AQN.j:end,AQN.j:end-1),-1)' * AQN.tau_j;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Compute Other parameters %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% Compute Omega
Omega = AQN.S(:,AQN.j:end)'*AQN.y_j*AQN.b_j' + AQN.S(:,AQN.j:end)'*AQN.Y(:,AQN.j:end-1) - triu(AQN.SY(AQN.j:end,AQN.j:end-1),0);

% Compute omega
omega = AQN.b_j/sqrt(AQN.rho_j);

% First: vectorize the upper part of the difference matrix
length = (size(AQN.Y,2)-AQN.j)*(size(AQN.Y,2)-AQN.j+1);
startPoint = zeros(length,1);
A_j_temp = AQN.mapToMatrix(startPoint,size(AQN.Y,2)-AQN.j);
vec_temp_1 = [];
for i = 1:size(A_j_temp,2)
    vec_temp_1 = [vec_temp_1 ; SHS_j(1:i,:)*A_j_temp(:,i) + AQN.b_j(i)*AQN.S(:,1:i)'*AQN.y_j];
end

% Second: vectorize the lower part of the diffrence matrix
func_temp_2 = A_j_temp' * SHS_j * A_j_temp + Omega'*A_j_temp + A_j_temp'*Omega - omega*omega';
vec_temp_2 = AQN.mapToVector(func_temp_2);

% Third: get the vectorization of the difference matrix
vec_temp = [vec_temp_1 ; vec_temp_2];

% Evaluate start norm of vec_temp
AQN.startValue = norm(vec_temp,inf);

% Set perturbation matrix
pertMatrix = AQN.findiffstep * eye(length);

% Initialize iteration number
iter = 0;

% Initialize flag information
flag = 1;

% Start writing into file
output_directory = 'output';
output_file_name = strcat(output_directory,'/Newton_',AQN.problem_name,num2str(AQN.iter_num),'.out');
file_id = fopen(output_file_name,'w');

% Write basic information
fprintf(file_id,'Dimension size of MA_j is %6d * %6d\n',size(A_j_temp,1),size(A_j_temp,2));
fprintf(file_id,'The condition number of Q matrix is %+14.6e\n',cond(SHS_j));
fprintf(file_id,'The accuracy for direct search is %+14.6e\n',AQN.DSacc);
fprintf(file_id,'\n');
fprintf(file_id,'\n');

% Write header
fprintf(file_id,'-------------------------\n');
fprintf(file_id,'\n');
fprintf(file_id,'%14s %14s %14s %14s %14s %14s %14s %14s %14s %14s %14s %14s\n','iteration','|f|-inf','|f|_2','FD-ref','min(eig(JJ))','min(eig(JJ_damp))','cond(JJ_damp)','damp-coeff','|d|-inf','dir_der','stepsize','Newton/Steep');

% Initialize damping coefficient
damping_coeff = AQN.lambda;

while ~AQN.checkAccuracy(vec_temp)
    
    if iter > AQN.maxIter
        flag = 2;
        break;
    end
    
    % count iterations
    iter = iter + 1;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%     % print information per 100 iterations
%     if round(iter/100) * 100 == iter
%         iter
%         norm(vec_temp,inf)
%         norm(vec_temp,inf)/AQN.startValue
%     end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%    

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%     % Compute Jacobian by finite difference
%     finDiffGrad = [];
%     for i = 1:length
%         pertPoint = startPoint + pertMatrix(:,i);
%         A_j_pert = AQN.mapToMatrix(pertPoint,size(AQN.Y,2)-AQN.j);
%         
%         vec_pert_1 = [];
%         for j = 1:size(A_j_pert,2)
%             vec_pert_1 = [vec_pert_1 ; SHS_j(1:j,:)*A_j_pert(:,j) + AQN.b_j(j)*AQN.S(:,1:j)'*AQN.y_j];
%         end
% 
%         func_pert_2 = A_j_pert' * SHS_j * A_j_pert + Omega'*A_j_pert + A_j_pert'*Omega - omega*omega';
%         vec_pert_2 = AQN.mapToVector(func_pert_2);
% 
%         vec_pert = [vec_pert_1 ; vec_pert_2];
%         finDiffGrad = [finDiffGrad (vec_pert-vec_temp)/AQN.findiffstep];
%     end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% 
    
    % Compute true Jacobian
    trueJaco = [];
    for i = 1:(size(AQN.Y,2)-AQN.j)
        trueJaco = [trueJaco ; [zeros(i,(size(AQN.Y,2)-AQN.j+1)*(i-1)) SHS_j(1:i,:) zeros(i,(size(AQN.Y,2)-AQN.j+1)*(size(AQN.Y,2)-AQN.j-i))]];
    end
    for i = 1:(size(AQN.Y,2)-AQN.j)
        for j = 1:i
            row1 = [zeros(1,(size(AQN.Y,2)-AQN.j+1)*(j-1)) (A_j_temp(:,i)'*SHS_j + Omega(:,i)') zeros(1,(size(AQN.Y,2)-AQN.j+1)*((size(AQN.Y,2)-AQN.j)-j))];
            row2 = [zeros(1,(size(AQN.Y,2)-AQN.j+1)*(i-1)) (A_j_temp(:,j)'*SHS_j + Omega(:,j)') zeros(1,(size(AQN.Y,2)-AQN.j+1)*((size(AQN.Y,2)-AQN.j)-i))];
            trueJaco = [trueJaco ; row1 + row2];
        end
    end  
    
    %%%%%%%%%%%%
    % Set finite difference gradient as true Jacobian
    finDiffGrad = trueJaco;
    %%%%%%%%%%%%
    
    % Increase damping coefficient if the damping matrix is not PD or is ill-conditioning
    while min(eig(trueJaco'*trueJaco + damping_coeff*eye(size(trueJaco,2)))) < AQN.lambda || cond(trueJaco'*trueJaco + damping_coeff*eye(size(trueJaco,2))) > 1e+15
        damping_coeff = damping_coeff * 10;
    end
    
    % compute search direction for Newton step
    d = -(trueJaco'*trueJaco + damping_coeff*eye(size(trueJaco,2)))\(trueJaco'*vec_temp);
    
    % set initial steplength
    alpha = 1;
    
    % set new point and update information needed (Newton Step)
    newPoint = startPoint + alpha*d;
    A_j_new = AQN.mapToMatrix(newPoint,size(AQN.Y,2)-AQN.j);
    
    % First: vectorize the upper part of the difference matrix
    vec_new_1 = [];
    for i = 1:size(A_j_new,2)
        vec_new_1 = [vec_new_1 ; SHS_j(1:i,:)*A_j_new(:,i) + AQN.b_j(i)*AQN.S(:,1:i)'*AQN.y_j];
    end

    % Second: vectorize the lower part of the difference matrix
    func_new_2 = A_j_new' * SHS_j * A_j_new + Omega'*A_j_new + A_j_new'*Omega - omega*omega';
    vec_new_2 = AQN.mapToVector(func_new_2);

    % Third: get the vectorization of the difference matrix
    vec_new = [vec_new_1 ; vec_new_2];
    
    % If sufficient decrease fails, do backtracking line search
    if norm(vec_new) >= norm(vec_temp)*(1-AQN.c1*alpha)
        
        while 1
            
            % Decrease the steplength and update information
            alpha = alpha/2;
            
            % set new point and update information needed (Newton Step)
            newPoint = startPoint + alpha*d;
            A_j_new = AQN.mapToMatrix(newPoint,size(AQN.Y,2)-AQN.j);
            
            % First: vectorize the upper part of the difference matrix
            vec_new_1 = [];
            for i = 1:size(A_j_new,2)
                vec_new_1 = [vec_new_1 ; SHS_j(1:i,:)*A_j_new(:,i) + AQN.b_j(i)*AQN.S(:,1:i)'*AQN.y_j];
            end

            % Second: vectorize the lower part of the difference matrix
            func_new_2 = A_j_new' * SHS_j * A_j_new + Omega'*A_j_new + A_j_new'*Omega - omega*omega';
            vec_new_2 = AQN.mapToVector(func_new_2);

            % Third: get the vectorization of the difference matrix
            vec_new = [vec_new_1 ; vec_new_2];
            
            % If sufficient decrease is satisfied
            if norm(vec_new) < norm(vec_temp)*(1-AQN.c1*alpha)
                break;
            end
            
            % If steplength reached minimum tolerance value, stop.
            if alpha < AQN.stepsize_limit
                norm(vec_new) - norm(vec_temp)
                norm(vec_temp)
                norm(vec_temp,Inf)
                fprintf('Newton method line search failed! \n');
                flag = 0;
                break;
            end
            
        end
        
    end
    
    % Set initial steplength of steepest descent
    alpha_steep = 1;
    
    % Set search direciton for steepest descent
    d_steep = -trueJaco'*vec_temp/norm(vec_temp);
    
    % set new point and update information needed (Steepest Descent Step)
    steepPoint = startPoint + alpha_steep * d_steep;
    A_j_steep = AQN.mapToMatrix(steepPoint,size(AQN.Y,2)-AQN.j);
    
    % First: vectorize the upper part of the difference matrix
    vec_steep_1 = [];
    for i = 1:size(A_j_steep,2)
        vec_steep_1 = [vec_steep_1 ; SHS_j(1:i,:)*A_j_steep(:,i) + AQN.b_j(i)*AQN.S(:,1:i)'*AQN.y_j];
    end

    % Second: vectorize the lower part of the difference matrix
    func_steep_2 = A_j_steep' * SHS_j * A_j_steep + Omega'*A_j_steep + A_j_steep'*Omega - omega*omega';
    vec_steep_2 = AQN.mapToVector(func_steep_2);

    % Third: get the vectorization of the difference matrix
    vec_steep = [vec_steep_1 ; vec_steep_2];
    
    % If sufficient decrease fails, do backtracking line search
    if norm(vec_steep) >= norm(vec_temp) - AQN.c1 * alpha_steep * (d_steep'*d_steep)
        
        while 1
            
            % Decrease the steplength and update information
            alpha_steep = alpha_steep/2;
            
            % set new point and update information needed (Steepest Descent Step)
            steepPoint = startPoint + alpha_steep * d_steep;
            A_j_steep = AQN.mapToMatrix(steepPoint,size(AQN.Y,2)-AQN.j);
            
            % First: vectorize the upper part of the difference matrix
            vec_steep_1 = [];
            for i = 1:size(A_j_steep,2)
                vec_steep_1 = [vec_steep_1 ; SHS_j(1:i,:)*A_j_steep(:,i) + AQN.b_j(i)*AQN.S(:,1:i)'*AQN.y_j];
            end

            % Second: vectorize the lower part of the difference matrix
            func_steep_2 = A_j_steep' * SHS_j * A_j_steep + Omega'*A_j_steep + A_j_steep'*Omega - omega*omega';
            vec_steep_2 = AQN.mapToVector(func_steep_2);

            % Third: get the vectorization of the difference matrix
            vec_steep = [vec_steep_1 ; vec_steep_2];
            
            % If sufficient decrease is satisfied
            if norm(vec_steep) < norm(vec_temp) - AQN.c1 * alpha_steep * (d_steep'*d_steep)
                break;
            end
            
            % If steplength reached minimum tolerance value, stop.
            if alpha_steep < AQN.stepsize_limit
                steepPoint = startPoint;
                A_j_steep = A_j_temp;
                vec_steep = vec_temp;
                alpha_steep = 0;
                break;
            end
            
        end
        
    end
    
    % Check which one performs better: steepest descent or Newton's method
    if norm(vec_new) < norm(vec_steep)
    
        % print information for each iteration
        fprintf(file_id,'%14d %+14.6e %+14.6e %+14.6e %+14.6e %+14.6e %+14.6e %+14.6e %+14.6e %14.6e %+14.6e %14s\n',iter,norm(vec_temp,inf),norm(vec_temp),max(max(abs(trueJaco-finDiffGrad))),min(eig(trueJaco'*trueJaco)),min(eig(trueJaco'*trueJaco + damping_coeff*eye(size(trueJaco,2)))),cond(trueJaco'*trueJaco + damping_coeff*eye(size(trueJaco,2))),damping_coeff,norm(d,inf),(vec_temp'*trueJaco*d)/norm(vec_temp),alpha,'Newton');

        % Update information
        startPoint = newPoint;
        A_j_temp = A_j_new;
        vec_temp = vec_new;
        
        % If steplength is less than 1, increase damping coefficient; otherwise, decrease damping coefficient
        if alpha < 1
            damping_coeff = damping_coeff * 2;
        else
            damping_coeff = max(damping_coeff/2,AQN.lambda);
        end
        
    else
        
        % print information for each iteration
        fprintf(file_id,'%14d %+14.6e %+14.6e %+14.6e %+14.6e %+14.6e %+14.6e %+14.6e %+14.6e %14.6e %+14.6e %14s\n',iter,norm(vec_temp,inf),norm(vec_temp),max(max(abs(trueJaco-finDiffGrad))),min(eig(trueJaco'*trueJaco)),min(eig(trueJaco'*trueJaco + damping_coeff*eye(size(trueJaco,2)))),cond(trueJaco'*trueJaco + damping_coeff*eye(size(trueJaco,2))),damping_coeff,norm(d_steep,inf),(trueJaco'*vec_temp/norm(vec_temp))'*(trueJaco'*vec_temp/norm(vec_temp)),alpha_steep,'steep');

        % Update information
        startPoint = steepPoint;
        A_j_temp = A_j_steep;
        vec_temp = vec_steep;
        
        % Increase damping coefficient
        damping_coeff = damping_coeff * 2;
        
    end
    
    % Terminate with Newton line search failed
    if flag == 0
        break;
    end
    
end

% Print two blank lines
fprintf(file_id,'\n');
fprintf(file_id,'\n');

% print termination information based on flag value
if flag == 1
    fprintf(file_id,'Success! \n');
elseif flag == 0
    fprintf(file_id,'Newton method line search failed! \n');
elseif flag == 2
    fprintf(file_id,'Iteration Limit reached! \n');
end

% Close file editting
fclose(file_id);

% Check whether accuracy reached
AQN.checkFlag = AQN.checkAccuracyNewton(vec_temp);

% Set result matrix MA_j
AQN.A_j = A_j_temp;

% Run overall unit tests
AQN.runUnitTestOld(5,SHS_j,vec_temp);


end
