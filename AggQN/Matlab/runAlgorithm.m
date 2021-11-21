function runAlgorithm(algorithm,history)

% Authors     : Albert Berahas, Frank E. Curtis, Baoyu Zhou
% Description : Runs algorithm to solve problem
% Input       : algorithm ~ algorithm name
%               history   ~ history length

% Set parameters
parameters.iteration_limit        = 1000;
parameters.cpu_limit              = 3600;
parameters.stationarity_tolerance = 1e-06;
parameters.damping_coefficient    = 1e-04;
parameters.c1                     = 1e-08;
parameters.c2                     = 5e-01;
parameters.stepsize_limit         = 1e-15;

% Set output file
file_id = create_output_file(algorithm,history);

% Read problem
x = [-2;2];

% Determine problem size
n = length(x);

% Determine scaling
parameters.scaling = 1;

% Initialize counters
counters.iteration = 1;
counters.objective = 0;
counters.gradient  = 0;

% Compute function and gradient
[f,counters] = evaluate_function(x,counters,parameters,0);
[g,counters] = evaluate_function(x,counters,parameters,1);

% Store initial gradient norm
norm_g_inf_initial = norm(g,inf);

% Set initial Hessian approximation
W_matr = speye(n);
W_func = @(x) x;
H_func = @(x) x;

% Initialize AggQN class
if strcmp(algorithm,'BFGS') == 1
  AQN = AggQN('denseInverseHessian',W_matr);
  AQN.setVerbosity(0);
elseif strcmp(algorithm,'LBFGS') == 1
  AQN = AggQN('limitedMemory',W_func,H_func,n,min(history,n));
  AQN.setVerbosity(0);
  AQN.setAggregation(false);
else
  AQN = AggQN('limitedMemory',W_func,H_func,n,min(history,n));
  AQN.setVerbosity(0);
  AQN.setAggregation(true);
end

% Print output header
print_header(n,algorithm,history,file_id);

% Initialize stepsize
alpha   = 1;
ls_flag = 'NA';

% Initialize messages
msg = '---';
AQN_msg = '---';

% Initialize timer
tic;

% Main loop
while 1

  % Save iteration information
  print_iteration(file_id,f,g,0,counters,alpha,ls_flag,AQN,AQN_msg);

  % Check for termination
  if norm(g,inf) <= parameters.stationarity_tolerance * max(norm_g_inf_initial,1.0)
    msg = 'Stationary point found';
    break;
  end
  if counters.iteration >= parameters.iteration_limit
    msg = 'Iteration limit reached';
    break;
  end
  if toc >= parameters.cpu_limit
    msg = 'CPU limit reached';
    break;
  end

  % Compute search direction
  d = AQN.computeInverseHessianProduct(-g);

  % Store current gradient
  g_stored = g;

  % Run line search
  [alpha,ls_flag,x,f,g,counters] = run_line_search(parameters,counters,x,d,f,g);

  % Check for line search failure
  if ls_flag(1) == 'F'
    msg = 'Line search failed';
    break;
  end

  % Compute (s,y) pair
  s = alpha*d;
  y = g - g_stored;

  % Damping
  if s'*y < parameters.damping_coefficient*(s'*s)
    theta = (1 - parameters.damping_coefficient)*(s'*s)/(s'*s - s'*y);
    y = theta*y + (1-theta)*s;
  end

  % Add pair
  AQN_msg = AQN.addPair(s,y);

  % Increment iteration counter
  counters.iteration = counters.iteration + 1;

end

% Print footer
print_footer(file_id,counters,f,g,msg,AQN);

% Close file
fclose(file_id);

end

function [v,counters] = evaluate_function(x,counters,parameters,option)

% Authors     : Albert Berahas, Frank E. Curtis, Baoyu Zhou
% Description : Evaluates problem function and increases counter
% Input       : x          ~ iterate
%               counters   ~ evaluation counters
%               parameters ~ parameters
%               option     ~ evaluation option
%                        0 ~ evaluate objective
%                        1 ~ evaluate gradient
% Output      : v          ~ return value
%               counters   ~ updated evaluation counters

% Check evaluation option
if option == 0

  % Increment counter
  counters.objective = counters.objective + 1;

  % Evaluate objective
  v = 100*(x(2)-x(1)^2)^2+(1-x(1))^2;

  % Scale objective value
  v = v/parameters.scaling;

else

  % Increment counter
  counters.gradient = counters.gradient + 1;

  % Evaluate gradient
  v = [-x(1)*400*(x(2)-x(1)^2)-2*(1-x(1));
             200*(x(2)-x(1)^2)            ];

  % Scale gradient value
  v = v/parameters.scaling;

end

end

function file_id = create_output_file(algorithm,history)

% Authors     : Albert Berahas, Frank E. Curtis, Baoyu Zhou
% Description : Creates output file
% Input       : algorithm ~ algorithm name
%               history   ~ history length
% Output      : file_id   ~ file ID for printing

% Set output_file_name
output_file_name = strcat('rosenbrock_',algorithm);
if strcmp(algorithm,'BFGS') ~= 1
  output_file_name = strcat(output_file_name,'_',num2str(history));
end
output_file_name = strcat(output_file_name,'.out');

% Set file id
file_id = fopen(output_file_name,'w');

end

function print_header(n,algorithm,history,file_id)

% Authors     : Albert Berahas, Frank E. Curtis, Baoyu Zhou
% Description : Prints output header
% Input       : n         ~ problem size
%               algorithm ~ algorithm name
%               history   ~ history length
%               file_id   ~ file ID for printing

% Print header
fprintf(file_id,'-------------------------\n');
fprintf(file_id,'Problem   : %13s\n','rosenbrock');
fprintf(file_id,'Dimension : %13d\n',n);
fprintf(file_id,'Algorithm : %13s\n',algorithm);
fprintf(file_id,'History   : %13d\n',history);
fprintf(file_id,'-------------------------\n');
fprintf(file_id,'\n');
fprintf(file_id,'%9s %8s %8s %14s %14s %14s %14s %8s %8s %14s %14s %10s %8s\n','iteration','# func.','# grad.','f','||g+A''lam||','||c||','stepsize','flag','# pairs','agg. index','agg. error','agg. count','agg. msg');

end

function print_iteration(file_id,f,g,c,counters,alpha,ls_flag,AQN,AQN_msg)

% Authors     : Albert Berahas, Frank E. Curtis, Baoyu Zhou
% Description : Prints output for an iteration
% Input       : file_id      ~ file ID for printing
%               f            ~ objective value
%               g            ~ gradient value
%               counters     ~ counters
%               alpha        ~ stepsize from last iteration
%               ls_flag      ~ line search flag from last iteration
%               AQN          ~ AggQN object

% Print iteration information
fprintf(file_id,'%9d %8d %8d %+14e %+14e %+14e %+14e %8s %8d %14d %+14e %10d %8s\n',counters.iteration,counters.objective,counters.gradient,f,norm(g,inf),norm(c,inf),alpha,ls_flag,AQN.numberOfPairs,AQN.aggregationIndex,AQN.aggregationError,AQN.aggregationCount,AQN_msg);

end

function print_footer(file_id,counters,f,g,msg,AQN)

% Authors     : Albert Berahas, Frank E. Curtis, Baoyu Zhou
% Description : Prints output for an iteration
% Input       : file_id  ~ file ID for printing
%               counters ~ counters from algorithm
%               f        ~ objective value
%               g        ~ gradient value
%               msg      ~ termination message

fprintf(file_id,'\n');
fprintf(file_id,'EXIT: %s\n',msg);
fprintf(file_id,'\n');
fprintf(file_id,'Objective            : %+14.6e\n',f);
fprintf(file_id,'Gradient norm        : %+14.6e\n',norm(g,inf));
fprintf(file_id,'Iterations           : %14d\n',counters.iteration);
fprintf(file_id,'Function evaluations : %14d\n',counters.objective);
fprintf(file_id,'Gradient evaluations : %14d\n',counters.gradient);
fprintf(file_id,'Aggregations         : %14d\n',AQN.aggregationCount);

end

function [alpha,flag,x_new,f_new,g_new,counters] = run_line_search(parameters,counters,x,d,f,g)

% Authors     : Albert Berahas, Frank E. Curtis, Baoyu Zhou
% Description : Runs line search
% Input       : parameters ~ parameters
%               counters   ~ counters
%               x          ~ iterate
%               d          ~ search direction
%               f          ~ objective value
%               g          ~ gradient value
% Output      : alpha      ~ stepsize
%               flag       ~ termination flag
%               x_new      ~ new iterate
%               f_new      ~ objective value at x_new
%               g_new      ~ gradient value at x_new
%               counters   ~ updated evaluation counters

% Initialize stepsize parameters
alpha    = 1;
alpha_lo = 0;
alpha_hi = 1000;

% Compute directional derivative
directional_derivative = g'*d;

% Line search loop
while 1

  % Set trial point
  x_new = x + alpha*d;

  % Check for small stepsize
  if alpha_hi - alpha_lo <= parameters.stepsize_limit

    % Try function and gradient evaluations
    try
      [f_new,counters] = evaluate_function(x_new,counters,parameters,0);
      [g_new,counters] = evaluate_function(x_new,counters,parameters,1);
    catch
      alpha = 0;
      flag  = 'F';
      x_new = x;
      f_new = f;
      g_new = g;
      return;
    end

    % Check for function nonincrease
    if f_new <= f
      flag = 'OK';
    else
      alpha = 0;
      flag  = 'F';
      x_new = x;
      f_new = f;
      g_new = g;
    end

    % Return
    return;

  end

  % Initialize boolean
  evaluation_success = true;

  % Try function evaluation
  try
    [f_new,counters] = evaluate_function(x_new,counters,parameters,0);
  catch
    evaluation_success = false;
  end

  % Check for evaluation success
  if evaluation_success

    % Set Armijo condition bool
    armijo_condition = (f_new <= f + parameters.c1*alpha*directional_derivative);

    % Check Armijo condition
    if armijo_condition

      % Try gradient evaluation
      try
        [g_new,counters] = evaluate_function(x_new,counters,parameters,1);
      catch
        evaluation_success = false;
      end

      % Check for evaluation success
      if evaluation_success

        % Check curvature condition
        if d'*g_new >= parameters.c2*directional_derivative

          % Success!
          flag = 'AW';

          % Return
          return;

        else

          % Update lower bound
          alpha_lo = alpha;

        end

      else

        % Update upper bound
        alpha_hi = alpha;

      end

    else

      % Update upper bound
      alpha_hi = alpha;

    end

  else

    % Update upper bound
    alpha_hi = alpha;

  end

  % Update stepsize
  if alpha_hi < inf
    alpha = 0.5*(alpha_lo + alpha_hi);
  else
    alpha = 2*alpha;
  end

end

end
