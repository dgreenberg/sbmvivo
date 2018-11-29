function X = randperm_jansimon(Arg1, Arg2, Arg3)  %#ok<INUSD,STOUT>
% Random permutation of array elements
% randperm_jansimon works in Inplace, Index and Derange mode:
%
% 1. INPLACE MODE: Y = randperm_jansimon(X, Dim)   - - - - - - - - - - - - - - - - - - - -
% The elements of the input array are randomly re-ordered along the specified
% dimension. In opposite to X(RANDPERM(LENGTH(X)) no temporary memory is used,
% such that randperm_jansimon is much more efficient for large arrays.
% INPUT:
%   X: Array of any size. Types: DOUBLE, SINGLE, CHAR, LOGICAL,
%      (U)INT64/32/16/8. The processed dimension must be shorter than 2^32
%      elements. To shuffle larger or complex arrays, use the Index mode.
%   Dim: Dimension to operate on.
%      Optional, default: [], operate on first non-singleton dimension.
% OUTPUT:
%   Y: Array of same size and type as X, but with shuffled elements.
%
% 2. INDEX MODE: Index = randperm_jansimon(N, 'index', NOut)   - - - - - - - - - - - - - -
% A vector of shuffled indices is created as by RANDPERM, but faster and using
% the smallest possible integer type to save memory. The number of output
% elements can be limited. This method works for cells, structs or complex
% arrays.
% INPUT:
%   N:      Numeric scalar >= 0.
%   String: 'index'.
%   NOut:   Optional, number of outputs, NOut <= N. Default: N.
% OUTPUT:
%   Index:  [1 x nOut] vector containing the shuffled elements of the vector
%           [1:N] vector. To limit the memory usage the smallest possible type
%           is used: UINT8/16/32 or INT64.
%
% 3. DERANGEMENT MODE: Index = randperm_jansimon(N, 'derange', NOut)   - - - - - - - - - -
% Equivalent to Index mode, but all elements Index[i] ~= i.
%
% CONTROL COMMANDS:  - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
% randperm_jansimon(S, 'seed'): S is either empty, scalar or a [1 x 4] DOUBLE vector to
%      set the random number  generator to a certain state. 4 integers between 0
%      and 2^32-1 are recommended for maximum entropy. If S is empty the default
%      seed is used.
% randperm_jansimon([], 'lock'): Lock Mex in the memory such that "clear all" does not
%      reset the random  number generator. If the Mex was compiled with
%      -DAUTOLOCK it is locked automatically.
% randperm_jansimon([], 'unlock'): Unlock the Mex file to allow a recompilation.
% randperm_jansimon(): Display compiler settings and status: Fast/exact integer, accept
%      cell/struct, current lock status, auto-lock mode.
%
% EXAMPLES:
%   randperm_jansimon(12345678901, 'seed');
%   R = randperm_jansimon(1:8)             %  [5, 8, 7, 4, 2, 6, 3, 1]
%   R = randperm_jansimon('abcdefg')       %  'gdfbcea'
%   R = randperm_jansimon([1:4; 5:8], 2)   %  [3, 4, 2, 1;  8, 6, 7, 5]
%   I = randperm_jansimon(8, 'index');     %  UINT8([3, 6, 5, 2, 4, 7, 8, 1])
% Choose 10 different rows from a 1000 x 100 matrix:
%   X = rand(1000, 100);       Y = X(randperm_jansimon(1000, 'index', 10), :);
% Operate on cells or complex arrays:
%   C = {9, 's', 1:5};         SC = C(randperm_jansimon(numel(C), 'index'));
%   M = rand(3) + i * rand(3); SM = M(:, randperm_jansimon(size(C, 2), 'index'))
%
% NOTES: randperm_jansimon(X) is about 50% to 85% faster than: Y = X(randperm(numel(X)).
%   It uses the Knuth-shuffle algorithm, also known as Fisher-Yates-shuffle.
%   The random numbers are created by the cute KISS algorithm of George
%   Marsaglia, which has a period of about 2^124 (~10^37).
%   More notes in randperm_jansimon.c.
%
% COMPILATION: See randperm_jansimon.c
%
% Tested: Matlab 6.5, 7.7, 7.8, WinXP, 32bit
%         Compiler: LCC2.4/3.8, BCC5.5, OWC1.8, MSVC2008
% Assumed Compatibility: higher Matlab versions, Mac, Linux, 64bit
% Author: Jan Simon, Heidelberg, (C) 2010-2011 j@n-simon.de
%
% See also RAND, RANDPERM.
% FEX:

% $JRev: R-m V:012 Sum:Q65/zXZqnvFX Date:07-Mar-2011 00:49:06 $
% $License: BSD (use/copy/change/redistribute on own risk, mention the author) $
% $UnitTest: uTest_Shuffle $
% $File: Tools\GLSets\Shuffle.m $

% History see randperm_jansimon.c

% Some Matlab algorithms for the Knuth-Shuffle: --------------------------------
% They are faster than Matlab's RANDPERM methods, but the MEX is recommended!

% n = numel(X);
% for i = 2:n      % Knuth shuffle in forward direction:
%    w    = ceil(rand * i);   % 1 <= w <= i
%    t    = X(w);
%    X(w) = X(i);
%    X(i) = t;
% end
 
% for i = n:-1:2   % Knuth shuffle in backward direction:
%    w    = ceil(rand * i);   % 1 <= w <= i
%    t    = X(w);
%    X(w) = X(i);
%    X(i) = t;
% end

% for i = 1:nOut   % Limit output:
%    w    = ceil(rand * (n - i + 1)) + (i - 1);   % i <= w <= n
%    t    = X(w);
%    X(w) = X(i);
%    X(i) = t;
% end
% X = X(1:nOut);

error(['JSimon:', mfilename, ':NoMex'], 'Need compiled mex file!');
