% Copyright (C) 2019 Albert Berahas, Frank E. Curtis, Baoyu Zhou
%
% All Rights Reserved.
%
% Authors: Albert Berahas, Frank E. Curtis, Baoyu Zhou
%
% Aggregated Quasi-Newton (AggQN) class definition
%
% Please cite:
%
%   A. Berahas, F. E. Curtis, and B. Zhou. "Limited-Memory BFGS with
%     Displacement Aggregation." arXiv, 1903.03471. 2019.
%
% Constructors:
%
%   AQN = AggQN('limitedMemory',initWv,initHv,n,m)
%         where storage mode 'limitedMemory' means
%               initWv is a handle of a function that computes
%                 matrix-vector products with the "initial" inverse Hessian,
%               initHv is a handle of a function that computes
%                 matrix-vector products with the "initial" Hessian,
%               n is the length of s and y vectors, and
%               m is the (limited memory) history length
%
%   AQN = AggQN({'denseInverseHessian','denseHessian'},M)
%         where storage mode 'denseInverseHessian' means
%               M is "initial" inverse Hessian
%         where storage mode 'denseHessian' means
%               M is "initial" Hessian
%
% Public methods:
%
%   msg = AQN.addPair(s,y)
%       Adds pair (s,y) to pairs defining (inverse) Hessian matrix.  If
%       storage mode is 'limitedMemory', then pair is added to matrices
%       (S,Y) and other quantities are updated based on other options.  If
%       storage mode is 'denseInverseHessian' or 'denseHessian', then the
%       inverse Hessian (resp., Hessian) matrix is updated using the pair
%       (s,y) as in standard BFGS.  Return value indicates result of adding
%       pair: 'Add' = pair added, 'Agg' = pair aggregated, or 'Swp' = pair
%       swapped.
%
%   e = AQN.aggregationCount
%       Returns counter of aggregations.  If storage mode is
%       'limitedMemory', then this value is nonnegative;
%       otherwise, it is -1.
%
%   e = AQN.aggregationError
%       Returns most recent error of aggregation.  If storage mode is
%       'limitedMemory' and an aggregation has occurred, then this value is
%       nonnegative;
%       otherwise, it is -1.
%
%   j = AQN.aggregationIndex
%       Returns most recent index of aggregated pair.  If storage mode is
%       'limitedMemory' and an aggregation has occured, then this value is
%       positive;
%       otherwise, it is -1.
%
%   Hv = AQN.computeHessianProduct(v)
%       Returns Hessian-vector product.  If storage mode is
%       'limitedMemory', then this is inefficient since the Hessian
%       approximate needs to be constructed (in compact form).  If storage
%       mode is 'denseInverseHessian', then this is inefficient since a
%       linear system needs to be solved.  If storage mode is
%       'denseHessian', then the matrix-vector product is computed
%       directly.
%
%   Wv = AQN.computeInverseHessianProduct(v)
%       Returns inverse-Hessian-vector product.  If storage mode is
%       'limitedMemory', then a two-loop recursion is employed.  If storage
%       mode is 'denseInverseHessian', then the matrix-vector product is
%       computed directly.  If storage mode is 'denseHessian', then this is
%       inefficient since a linear system needs to be solved.
%
%   H = AQN.Hessian
%       Returns Hessian approximation.  If storage mode is 'limitedMemory'
%       or 'denseInverseHessian', then this is inefficient and should be
%       avoided; it would likely be better to compute
%       inverse-Hessian-vector products than call this function.
%
%   W = AQN.inverseHessian
%       Returns inverse Hessian approximation.  If storage mode is
%       'limitedMemory' or 'denseHessian', then this is inefficient and
%       should be avoided; it would likely be better to compute
%       inverse-Hessian-vector products than call this function.
%
%   p = AQN.numberOfPairs
%       Returns number of pairs currently stored.
%
%   AQN.printData
%       Prints all members of the class.
%
%   AQN.setAccuracyTolerance(value)
%       Sets accuracy tolerance for aggregation to value.
%
%   AQN.setAggregation(b)
%       Sets aggregation to b=true (default) or b=false.  Storage mode must
%       be 'limitedMemory' or an error is thrown.  If aggregation is true,
%       then gradient displacements may be modified.  In non-aggregation
%       mode, gradient displacements are not modified, as in traditional
%       L-BFGS.
%
%   AQN.setDebug(b)
%       Sets debug option to b=true or b=false (default).  In debug mode,
%       unit tests are performed to check the correctness and accuracy of
%       intermediate operations.  If the verbosity level is equal to 0 when
%       debug is set to true, then the verbosity level is increased to 1.
%       In non-debug mode, the unit tests are not run.
%
%   AQN.setOnlySwap(b)
%       Sets onlySwap to b=true or b=false (default).  If aggregation is
%       true and onlySwap is false, then aggregation is performed.  If
%       aggregation is true and onlySwap is true, then the same strategy is
%       used to identify which pairs to remove (rather than always the
%       oldest, as in standard L-BFGS), but gradient displacements are not
%       modified.
%
%   AQN.tryNewton(b)
%       Sets tryNewton to b=true or b=false (default).  If aggregation is
%       true and tryNewton is true, then a Newton iteration is performed if
%       the error in the aggregation is sufficiently large.  This can be
%       inefficient and should only be used for testing purposes!
%
%   AQN.setVerbosity(level)
%       Sets verbosity to level, where level must be a nonnegative integer.
%       - verbosity=0 means that nothing is printed (default).
%       - verbosity=1 means that warnings and short messages are printed.
%       - verbosity=2 means that all data is printed after pair is added.

% AggQN class
classdef AggQN < handle
  
  % Properties (private access)
  properties (SetAccess = private, GetAccess = private)
    
    %%%%%%%%%%%
    % Options %
    %%%%%%%%%%%
    storage_mode = 'limitedMemory' % storage mode
                                   %  'limitedMemory'       = (S,Y) is stored; H and W are empty
                                   %  'denseInverseHessian' = W     is stored; (S,Y) and H are empty
                                   %  'denseHessian'        = H     is stored; (S,Y) and W are empty
    aggregate    = true            % aggregation mode
    debug        = false           % debug mode
    onlySwap     = false           % only swap in aggregation
    tryNewton    = false           % try Newton mode
    verbosity    = 0               % verbosity level
                                   %  0 = prints nothing
                                   %  1 = prints warnings and short messages
                                   %  2 = also prints data after pair added
    
    %%%%%%%%%%%%
    % Counters %
    %%%%%%%%%%%%
    count = 0

    %%%%%%%%%%%%%%%%%%%%%%%%%%%%
    % (Inverse) Hessian values %
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%
    m = inf % history length
    n = -1  % Vector length / matrix size
    S = []  % Iterate displacements
    Y = []  % Gradient displacements
    W = []  % Inverse Hessian approximation
    H = []  % Hessian approximation
    
    %%%%%%%%%%%%%%%%%%%%
    % Function handles %
    %%%%%%%%%%%%%%%%%%%%
    initWv = [] % "Initial" inverse Hessian product function handle
    initHv = [] % "Initial" Hessian product function handle
    
    %%%%%%%%%%%%%%%%%%%%
    % Auxiliary values %
    %%%%%%%%%%%%%%%%%%%%
    rho = [] % Reciprocal of displacement products
    SY  = [] % S'*Y
             %   where S = [s_earliest ... s_latest]
             %     and Y = [y_earliest ... y_latest]
    R   = [] % Compact form matrix; see Byrd, Nocedal, and Schnabel (1994)
    D   = [] % Compact form matrix; see Byrd, Nocedal, and Schnabel (1994)
    L   = [] % Compact form matrix; see Byrd, Nocedal, and Schnabel (1994)
    HS  = [] % initHv(S)
    SHS = [] % S'*initHv(S)
    
    %%%%%%%%%%%%%%%%%%%%%%
    % Aggregation values %
    %%%%%%%%%%%%%%%%%%%%%%
    j = -1
    s_j
    y_j
    rho_j
    tau_j
    HS_j
    SHS_j
    invM_j
    Omega_j
    omega_j
    rhs_j
    b_j
    MA_j
    A_j
    error = 0
    
    %%%%%%%%%%%%%%
    % TOLERANCES %
    %%%%%%%%%%%%%%
    chol_pert_init   = 1e-15
    proj_tol         = 1e-08
    proj_tol_loose   = 1e-04
    proj_tol_beta    = 1e-08
    curv_tol         = 1e-04
    acc_tol          = 1e-06
    acc_tol_newton   = 1e-06
    max_iter_newton  = 1e+01
    c1_newton        = 1e-08
    alpha_min_newton = 1e-12
    reg_newton       = 1e-14
    
  end
  
  % Methods (private access)
  methods (Access = private)
    
    % Add data, storage mode 'limitedMemory'
    addDataLimitedMemory(AQN,s,y,reason)
    
    % Add pair, storage mode 'denseHessian'
    msg = addPairDenseHessian(AQN,s,y)
    
    % Add pair, storage mode 'limitedMemory'
    msg = addPairLimitedMemory(AQN,s,y)
    
    % Add pair, storage mode 'denseInverseHessian'
    msg = addPairDenseInverseHessian(AQN,s,y)
    
    % Cholesky factorization, perturbed when needed
    [L,perturbation] = choleskyPerturb(AQN,M)
    
    % Compute error in solving aggregation equations
    [errorMaxAbs,errorFroSqr,v1,v2] = computeAggregationValuesError(AQN,A_j)

    % Compute remainder of aggregation values using exact method
    computeAggregationValuesExact(AQN)
    
    % Initialize aggregation values
    computeAggregationValuesInitial(AQN)

    % Compute aggregation values using Newton's method
    computeAggregationValuesNewton(AQN)
    
    % Compute inner Hessian product for aggregation
    Hjv = computeInnerHessianProduct(AQN,v)
    
    % Compute inner inverse Hessian product for aggregation
    Wjv = computeInnerInverseHessianProduct(AQN,v)
    
    % Constructor, storage mode 'denseHessian'
    constructorDenseHessian(AQN,H)
    
    % Constructor, storage mode 'limitedMemory'
    constructorLimitedMemory(AQN,initWv,initHv,n,m)
    
    % Constructor, storage mode 'denseInverseHessian'
    constructorDenseInverseHessian(AQN,W)
    
    % Delete data, storage mode 'limitedMemory'
    deleteDataLimitedMemory(AQN,s,y,reason)
    
    % Run unit test
    runUnitTest(AQN,number,level,phi,phi_comb,N)
            
  end
  
  % Methods (public access)
  methods (Access = public)
    
    % Constructor
    function AQN = AggQN(varargin)
      
      % Call specific constructor
      if nargin == 5 && strcmp(varargin{1},'limitedMemory') == 1
        AQN.storage_mode = 'limitedMemory';
        AQN.constructorLimitedMemory(varargin{2},varargin{3},varargin{4},varargin{5});
      elseif nargin == 2 && strcmp(varargin{1},'denseInverseHessian') == 1
        AQN.storage_mode = 'denseInverseHessian';
        AQN.constructorDenseInverseHessian(varargin{2});
      elseif nargin == 2 && strcmp(varargin{1},'denseHessian') == 1
        AQN.storage_mode = 'denseHessian';
        AQN.constructorDenseHessian(varargin{2});
      else
        msg = 'AggQN: Incorrect number of inputs for given storage mode';
        msg = [msg newline AQN.usageMessage];
        error(msg);
      end
      
    end
    
    % Add pair
    msg = addPair(AQN,s,y)
    
    % Aggregation count
    function c = aggregationCount(AQN)
      
      % Set aggregation count
      c = AQN.count;
      
    end
    
    % Aggregation error
    function e = aggregationError(AQN)
      
      % Set aggregation error
      e = AQN.error;
      
    end
        
    % Aggregated index
    function j = aggregationIndex(AQN)
      
      % Set as index of most recently aggregated pair
      j = AQN.j;
      
    end
    
    % Compute Hessian-vector product (H*v)
    Hv = computeHessianProduct(AQN,v)
    
    % Compute inverse-Hessian-vector product (W*v)
    Wv = computeInverseHessianProduct(AQN,v)
    
    % Hessian
    H = hessian(AQN)
    
    % Inverse Hessian
    W = inverseHessian(AQN)
    
    % Number of pairs
    function p = numberOfPairs(AQN)
      
      % Return number of pairs if storage mode is 'limitedMemory'
      if strcmp(AQN.storage_mode,'limitedMemory') == 1
        p = size(AQN.S,2);
      else
        p = -1;
      end
      
    end
    
    % Print data (for debugging)
    printData(AQN)
    
    % Sets accuracy tolerance
    setAccuracyTolerance(AQN,value)
    
    % Sets aggregation
    setAggregation(AQN,value)
    
    % Sets debug
    setDebug(AQN,value)
    
    % Sets onlySwap
    setOnlyRemove(AQN,value)
    
    % Sets tryNewton
    setTryNewton(AQN,value)
    
    % Sets verbosity
    setVerbosity(AQN,value)
    
  end
  
  % Methods (static)
  methods (Static)
    
    % Usage message
    msg = usageMessage
    
  end
  
end