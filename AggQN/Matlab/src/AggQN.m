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
%   A. Berahas, F. E. Curtis, and B. Zhou, "Limited-Memory BFGS with
%   Displacement Aggregation: Foundations and Superlinear Convergence Rate
%   Guarantees," Lehigh ISE/COR@L Technical Report, 2019.
%
% Constructors:
%
%   AQN = AggQN('SY',initWv,initHv,n,m)
%         where storage mode 'SY'
%               means initWv is a handle of a function that computes
%               matrix-vector products with the "initial" inverse Hessian
%               and initHv is a handle of a function that computes
%               matrix-vector products with the "initial" Hessian
%               and the history m must be a positive integer or inf
%
%   AQN = AggQN({'W','H'},M)
%         where storage mode 'W'
%               means M is "initial" inverse Hessian
%         where storage mode 'H'
%               means M is "initial" Hessian
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
%       aggregated, 'Swp' = pair swapped.
%
%   j = AQN.aggregatedIndex
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
%   AQN.setVerbosity(level)
%       Sets verbosity to level, where level must be a nonnegative integer.
%       Verbosity=0 means that nothing is printed.
%       Verbosity=1 means that warnings and short messages are printed.
%       Verbosity=2 means that all data is printed after pair is added.

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
    aggregate    = true  % aggregation
                         %  true  = aggregation on
                         %  false = aggregation off
    debug        = true  % debug mode
    verbosity    = 0     % verbosity level
                         %  0 = prints nothing
                         %  1 = prints warnings and short messages
                         %  2 = also prints data after pair added
    precondition = 1     % precondition level
                         %  0 = not use preconditioner
                         %  1 = use preconditioner
    
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
    rho  = []  % Reciprocal of displacement products
    SY   = []  % S'*Y
               %   where S = [s_earliest ... s_latest]
               %     and Y = [y_earliest ... y_latest]
    L_SS = []  % Cholesky factor (lower triangular) of S'*S
               %   where S = [s_latest ... s_earliest]
    L_SHS = [] % Cholesky factor of S'*H*S
               %   where H = "initial" Hessian
    
    %%%%%%%%%%%%%%%%%%%%%%
    % Aggregation values %
    %%%%%%%%%%%%%%%%%%%%%%
    j = -1
    s_j
    y_j
    rho_j
    tau_j
    MA_j
    A_j
    b_j
    HS_j
    HQ_j
    C
    Q
    sum_diff = []
    
    %%%%%%%%%%%%%%
    % TOLERANCES %
    %%%%%%%%%%%%%%
    chol_pert_init = 1e-15
    cond_tol_1     = 1e+14
    cond_tol_2     = 1e+14
    cond_tol_3     = 1e+06
    cond_tol_beta  = 1e+14
    %lin_ind_tol    = 1e-08
    parallel_tol   = 1e-15
    
  end
  
  % Methods (private access)
  methods (Access = private)
    
    % Add data, storage mode 'SY'
    addDataSY(AQN,s,y)
    
    % Add pair, storage mode 'SY'
    msg = addPairSY(AQN,s,y)

    % Add pair, storage mode 'W'
    msg = addPairW(AQN,s,y)

    % Add pair, storage mode 'H'
    msg = addPairH(AQN,s,y)
    
    % Cholesky factorization, perturbed when needed
    [L,perturbation] = choleskyPerturb(AQN,M)
    
    % Compute aggregation values 'A' and 'b'
    computeAggregationValues(AQN)
    
    % Compute aggregation values 'A' and 'b' without preconditioners
    computeAggregationValuesOld(AQN)
    
    % Delete data, storage mode 'SY'
    deleteDataSY(AQN,s,y)
    
    % Run unit test
    runUnitTest(AQN,number,invQ,rhs,level,phi,phi_comb,N)
    
    % Run unit test
    runUnitTestOld(AQN,number,invQ,rhs,level,phi,phi_comb,N)
        
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
    
    % Aggregated index
    function j = aggregatedIndex(AQN)
      
      % Set as index of most recently aggregated pair
      j = AQN.j;
      
    end
    
    % Get S and Y pairs information
    function [S Y] = aggregatedSY(AQN)
        
        % Set as current s,y pairs
        S = AQN.S;
        Y = AQN.Y;
        
    end
    
    % Get s_j and y_j after rotated
    function [s_j y_j] = rotatedpair(AQN)
        
        % Set as the pair after rotated
        s_j = AQN.s_j;
        y_j = AQN.y_j;
        
    end
    
    function [sum_diff] = get_diff(AQN)
        sum_diff = AQN.sum_diff;      
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
      
      % Set as size of S
      p = size(AQN.S,2);
      
    end
    
    % Print data (for debugging)
    printData(AQN)
    
    % Sets adaptivity
    setAdaptivity(AQN,value)
    
    % Sets aggregation
    setAggregation(AQN,value)
    
    % Sets debug
    setDebug(AQN,value)
    
    % Sets verbosity
    setVerbosity(AQN,value)
    
    % Sets preconditioner level
    setPreconditioner(AQN,value)
    
  end
  
  % Methods (static)
  methods (Static)
    
    % Usage message
    msg = usageMessage
    
  end
  
end