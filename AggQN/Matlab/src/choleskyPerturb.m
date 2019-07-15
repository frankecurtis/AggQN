% Copyright (C) 2019 Albert Berahas, Frank E. Curtis, Baoyu Zhou
%
% All Rights Reserved.
%
% Authors: Albert Berahas, Frank E. Curtis, Baoyu Zhou
%
% Method definition for AggQN class

function [L,perturbation] = choleskyPerturb(AQN,M)

% Compute factorization
flag = 1;
perturbation = 0;
while flag ~= 0
  [L,flag] = chol(M + perturbation*eye(size(M,2)),'lower');
  if flag ~= 0 && perturbation == 0
    perturbation = AQN.chol_pert_init;
  else
    perturbation = perturbation * 10;
  end
end
