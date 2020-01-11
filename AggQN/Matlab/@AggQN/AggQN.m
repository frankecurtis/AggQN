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
%   AQN = AggQN('SY',initWv,initHv,n,m)
%         where storage mode 'SY' means
%               initWv is a handle of a function that computes
%                 matrix-vector products with the "initial" inverse Hessian,
%               initHv is a handle of a function that computes
%                 matrix-vector products with the "initial" Hessian,
%               n is the length of s and y vectors, and
%               m is the history length (for limited memory methods)
%
%   AQN = AggQN({'W','H'},M)
%         where storage mode 'W' means
%               M is "initial" inverse Hessian
%         where storage mode 'H' means
%               M is "initial" Hessian
%
% Public methods:
%
%   msg = AQN.addPair(s,y)
%       Adds pair (s,y) to pairs defining (inverse) Hessian matrix.  If
%       storage mode is 'SY', then pair is added to matrices (S,Y), then
%       matrices are updated based on adaptive and aggregate options.  If
%       storage mode is 'W' ('H'), then inverse Hessian (Hessian) matrix is
%       updated using the pair (s,y) as in standard BFGS.  Return value
%       indicates result of adding pair: 'Add' = pair added, 'Agg' = pair
%       aggregated, or 'Swp' = pair swapped.
%
%   e = AQN.aggregationCount
%       Returns counter of aggregations.  If storage mode is 'SY', then
%       this value is nonnegative;
%       otherwise, it is -1.
%
%   e = AQN.aggregationError
%       Returns most recent error of aggregation.  If storage mode is
%       'SY' and an aggregation has occurred, then this value is nonnegative;
%       otherwise, it is -1.
%
%   j = AQN.aggregationIndex
%       Returns most recent index of aggregated pair.  If storage mode is
%       'SY' and an aggregation has occurred, then this value is positive;
%       otherwise, it is -1.
%
%   Hv = AQN.computeHessianProduct(v)
%       Returns Hessian-vector product.  If storage mode is 'SY', then this
%       is inefficient since the Hessian matrix is constructed.  If storage
%       mode is 'W', then this is inefficient since a linear system is
%       solved.  If storage mode is 'H', then matrix-vector product is
%       computed.
%
%   Wv = AQN.computeInverseHessianProduct(v)
%       Returns inverse-Hessian-vector product.  If storage mode is 'SY',
%       then two-loop recursion is applied.  If storage mode is 'W', then
%       matrix-vector product is computed.  If storage mode is 'H', then
%       this is inefficient since a linear system is solved.
%
%   H = AQN.Hessian
%       Returns Hessian approximation.  If storage mode is 'SY' or 'W',
%       then this is inefficient to construct and should be avoided; it
%       would likely be better to compute inverse Hessian products.
%
%   W = AQN.inverseHessian
%       Returns inverse Hessian approximation.  If storage mode is 'SY',
%       then this is inefficient to construct and should be avoided; it
%       would likely be better to compute inverse Hessian products.
%       Similarly, if storage mode is 'H', then this is inefficient to
%       construct and should be avoided; it would likely be better to
%       compute Hessian products.
%
%   p = AQN.numberOfPairs
%       Returns number of pairs currently stored.
%
%   AQN.printData
%       Prints all members of the class.
%
%   AQN.setAdaptivity(b)
%       Sets adaptivity to b=true (default) or b=false.  Storage mode must
%       be 'SY' or an error is thrown.  In adaptive mode, the # of pairs
%       may increase or decrease.  In non-adaptive mode, the # of pairs
%       increases to history m, then never increases or decreases.
%
%   AQN.setAggregation(b)
%       Sets aggregation to be b=true (default) or b=false.  Storage mode
%       must be 'SY' or an error is thrown.  In aggregation mode, gradient
%       displacements may be modified.  In non-aggregation mode, gradient
%       displacements are not modifed, as in traditional L-BFGS.
%
%   AQN.setDebug(b)
%       Sets debug option to be b=true (default) or b=false.  In debug
%       mode, unit tests are performed to check the correctness and
%       accuracy of intermediate operations.  If the verbosity level is
%       equal to 0 when debug is set to true, then the verbosity level is
%       increased to 1.  In non-debug mode, the unit tests are not run.
%
%   AQN.setPrecondition(b)
%       Sets precondition to be b=true (default) or b=false.  Storage mode
%       must be 'SY' or an error is thrown.  Only affects aggregation mode.
%
%   AQN.setVerbosity(level)
%       Sets verbosity to level, where level must be a nonnegative integer.
%       - verbosity=0 means that nothing is printed.
%       - verbosity=1 means that warnings and short messages are printed.
%       - verbosity=2 means that all data is printed after pair is added.

% AggQN class
classdef AggQN < handle
  
  % Properties (private access)
  properties (SetAccess = private, GetAccess = private)
    
    %%%%%%%%%%%
    % Options %
    %%%%%%%%%%%
    storage_mode = 'SY'  % storage mode
                         %  'SY' = (S,Y) is stored; H and W are empty
                         %  'W'  = W     is stored; (S,Y) and H are empty
                         %  'H'  = H     is stored; (S,Y) and W are empty
    adaptive     = true  % adaptivity
                         %  true  = # pairs may increase or decrease
                         %  false = # pairs increases up to,
                         %          then remains at m
    aggregate    = true  % aggregation mode
    debug        = true  % debug mode
    precondition = false % preconditioning mode
    tryNewton    = false % try Newton mode
    verbosity    = 1     % verbosity level
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
    rho   = [] % Reciprocal of displacement products
    SY    = [] % S'*Y
               %   where S = [s_earliest ... s_latest]
               %     and Y = [y_earliest ... y_latest]
    
    %%%%%%%%%%%%%%%%%%%%%%
    % Aggregation values %
    %%%%%%%%%%%%%%%%%%%%%%
    j = -1
    s_j
    y_j
    rho_j
    tau_j
    H_j
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
    proj_tol_loose   = 1e-01
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
    
    % Add data, storage mode 'SY'
    addDataSY(AQN,s,y,reason)
    
    % Add pair, storage mode 'H'
    msg = addPairH(AQN,s,y)
    
    % Add pair, storage mode 'SY'
    msg = addPairSY(AQN,s,y)
    
    % Add pair, storage mode 'W'
    msg = addPairW(AQN,s,y)
    
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
    
    % Constructor, storage mode 'H'
    constructorH(AQN,H)
    
    % Constructor, storage mode 'SY'
    constructorSY(AQN,initWv,initHv,n,m)
    
    % Constructor, storage mode 'W'
    constructorW(AQN,W)
    
    % Delete data, storage mode 'SY'
    deleteDataSY(AQN,s,y,reason)
    
    % Run unit test
    runUnitTest(AQN,number,level,phi,phi_comb,N)
            
  end
  
  % Methods (public access)
  methods (Access = public)
    
    % Constructor
    function AQN = AggQN(varargin)
      
      % Call specific constructor
      if nargin == 5 && strcmp(varargin{1},'SY') == 1
        AQN.storage_mode = 'SY';
        AQN.constructorSY(varargin{2},varargin{3},varargin{4},varargin{5});
      elseif nargin == 2 && strcmp(varargin{1},'W') == 1
        AQN.storage_mode = 'W';
        AQN.constructorW(varargin{2});
      elseif nargin == 2 && strcmp(varargin{1},'H') == 1
        AQN.storage_mode = 'H';
        AQN.constructorH(varargin{2});
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
      
      % Return number of pairs if storage mode is 'SY'
      if strcmp(AQN.storage_mode,'SY') == 1
        p = size(AQN.S,2);
      else
        p = -1;
      end
      
    end
    
    % Print data (for debugging)
    printData(AQN)
    
    % Sets accuracy tolerance
    setAccuracyTolerance(AQN,value)
    
    % Sets adaptivity
    setAdaptivity(AQN,value)
    
    % Sets aggregation
    setAggregation(AQN,value)
    
    % Sets debug
    setDebug(AQN,value)
    
    % Sets preconditioner level
    setPrecondition(AQN,value)
    
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