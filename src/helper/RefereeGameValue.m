%%  REFEREEGAMEVALUE    Computes the maximum value of a quantum game
%
%   Documentation coming soon (maybe).

function rgval = RefereeGameValue(p,V,varargin)

    % set optional argument defaults: LVL=1
    [lvl] = opt_args({ 1 },varargin{:});

    % Get some basic values.
    [ma,mb] = size(p);
    oa = size(V,1);
    ob = size(V,2);
    rdim = size(V,5);

    % Now run the SDP that computes the maximum.
    cvx_begin quiet
    
        variable q(oa,ob,ma,mb,rdim,rdim) complex;
            
        maximize real(sum(sum(p.*squeeze(sum(sum(sum(sum(V.*q,6),5),1),2)))))
            
        subject to
          NPAHierarchyRef(q,lvl) == 1;
    cvx_end
    
    save q
    assignin('base','q',q);

    rgval = cvx_optval;

    % Deal with error messages.
    if(strcmpi(cvx_status,'Inaccurate/Solved'))
        warning('RefereeGameValue:NumericalProblems','Minor numerical problems encountered by CVX. Consider adjusting the tolerance level TOL and re-running the script.');
    elseif(strcmpi(cvx_status,'Inaccurate/Infeasible'))
        warning('RefereeGameValue:NumericalProblems','Minor numerical problems encountered by CVX. Consider adjusting the tolerance level TOL and re-running the script.');
    elseif(strcmpi(cvx_status,'Unbounded') || strcmpi(cvx_status,'Inaccurate/Unbounded') || strcmpi(cvx_status,'Failed'))
        error('RefereeGameValue:NumericalProblems',strcat('Numerical problems encountered (CVX status: ',cvx_status,'). Please try adjusting the tolerance level TOL.'));
    end
end