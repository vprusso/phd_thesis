function [ ok ] = IsClose( v1, v2 )
%--------------------------------------------------------------------------
% Syntax:       ok = IsClose( v1, v2 );
%
% Inputs:       v1 - first value
%               v2 - second value
%
% Outputs:      ok - a Boolean variable which checks if v1 and v2 are close
%                    enough to each other in terms of eps = 1e-8
% 
% Description:  Determines if two values close to each other.
%
%               Notes:
%                   * IsClose(NaN,NaN) == 0
%                   * IsClose(Inf,Inf) == 0
%                   * The function compares v1/v2-1, v2/v1-1 and v1-v2 to 
%                     "close enough"
%                   * If for all matrix elements at least one of the above 
%                     is close enough, the matrices are considered close 
%                     enough.
%                   * A bit heavy in implementation, so don't use in tight 
%                     inner loops/searches 
%
% Requires:     IsEql.m
%
% References:   [1] QLib : http://www.tau.ac.il/~quantum/qlib/qlib.html
%
% Author:       Vincent Russo (Adapted from Adapted from QLib_v_1_0 [1])
%
% Date:         June 9, 2014
%--------------------------------------------------------------------------

ok = 0;
eps = 1e-8;
if ~isequal(size(v1),size(v2))
    return; % Unequal dimensions
end

warning off MATLAB:divideByZero
the_ratio_1 = abs(v1./v2 - 1);
the_ratio_2 = abs(v2./v1 - 1);
warning on  MATLAB:divideByZero

the_diff  = abs(v1 - v2);

the_ratio = min(the_ratio_1,the_ratio_2);

delta = min(the_ratio, the_diff);

problem_places = [find(delta > eps) find(isnan(delta))];

ok = isempty(problem_places);


end

