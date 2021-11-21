function testAggBFGS

% Authors     : Albert Berahas, Frank E. Curtis, Baoyu Zhou
% Description : Runs algorithms

% Set algorithms
algorithms = {"BFGS", "LBFGS", "AggBFGS"};
histories  = [3;5];

% Loop over algorithms
for j = 1:length(algorithms)

  % Print solve message
  fprintf('Running with %7s... ',algorithms{j});

  % Check algorithm name
  if strcmp(algorithms{j},'BFGS') == 1

    % Run algorithm
    runAlgorithm(algorithms{j},inf);

  else

    % Loop over AggBFGS history lengths
    for k = 1:length(histories)

      % Run algorithm
      runAlgorithm(algorithms{j},histories(k));

    end

  end

  % Print done message
  fprintf('done\n');

end
