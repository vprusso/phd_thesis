function v = ip( A,B )
%--------------------------------------------------------------------------
% Syntax:       v = ip( A, B );
%
% Inputs:       A,B - pure states
%
% Outputs:      v - the inner product of vectors A and B
%
% Description:  Simple inner product of two operators
%
% Requires:     Nothing
%--------------------------------------------------------------------------

v = trace(A'*B);

end