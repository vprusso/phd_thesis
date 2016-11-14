%--------------------------------------------------------------------------
% 	NPAR_Classical_Value: 	Computes the classical value of the 
%							monogamy-of-entanglement game.
%
%	This function has 2 required arguments:
%		R: a cell array consisting of the bases for the referee.	
%		REPS: the number of parallel repetitions performed.
%
% 	MAX_CVAL = NPAR_CLASSICAL_VALUE(R,REPS) is the classical value for the
% 	monogamy-of-entanglement game.
%		           
% Requires:   CVX [3], QETLAB
%
% References: [1] "A convergent hierarch of semidefinite programs
%                  characterizing the set of quantum correlations" - M.
%                  Navascues, S. Pironio, A. Acin. 
%
%             [2] "A monogamy-of-entanglement game with applications to
%                  device-independent quantum cryptography - M. Tomamichel,
%                  S. Fehr, J. Kaniewski, S. Wehner.
%
%             [3] CVX - (http://cvxr.com/cvx/)
%--------------------------------------------------------------------------
function [max_cval,rho]  = MonogamyGameUnentangledValue( R, varargin )

% set optional argument defaults: REPS = 1
[reps] = opt_args({ 1 },varargin{:});

% Get some basic values and make sure that the input vectors are column
% vectors.
num_inputs = length(R);
num_outputs = length(R{1});
[xdim,ydim] = size(R{1}{1});

% Now tensor things together if we are doing more than 1 repetition
if(reps > 1)
    i_ind = zeros(1,reps);
    j_ind = zeros(1,reps);

    for i = 1:num_inputs^reps
        for j = 1:num_outputs^reps
            for l = reps:-1:1
                to_tensor{l} = R{i_ind(l)+1}{j_ind(l)+1};
            end
            newR{i}{j} = Tensor(to_tensor);

            j_ind = update_odometer(j_ind,num_outputs*ones(1,reps));
        end
        i_ind = update_odometer(i_ind,num_inputs*ones(1,reps));
    end
    R = newR;

    % Recalculate.
    num_inputs = length(R);
    num_outputs = length(R{1});
    [xdim,ydim] = size(R{1}{1});
end

max_cval = 0; 
R_win_sum = 0; 

if num_inputs == 2
    for i = 1:num_outputs
        for j = 1:num_outputs
            R_win_sum = R{1}{i} + R{2}{j};   
            cvx_begin sdp quiet 
                %cvx_precision best; 
                %#ok<*VUNUS>    % suppress MATLAB warnings for equality checks in CVX
                %#ok<*EQEFF>    % suppress MATLAB warnings for inequality checks in CVX            
                variable rho(xdim,ydim) hermitian

                maximize trace(R_win_sum'*rho)
                subject to
                    trace(rho) == 1;
                    rho >= 0;        
            cvx_end

            cval = 1/num_inputs * cvx_optval;            

            if cval > max_cval
                max_cval = cval;            
            end
        end
    end
end
if num_inputs == 3    
    
    for i = 1:num_outputs
        for j = 1:num_outputs
            for k = 1:num_outputs
                
                R_win_sum = R{1}{i} + R{2}{j} + R{3}{k};

                cvx_begin sdp quiet 
                    %cvx_precision best; 
                    %#ok<*VUNUS>    % suppress MATLAB warnings for equality checks in CVX
                    %#ok<*EQEFF>    % suppress MATLAB warnings for inequality checks in CVX            
                    variable rho(xdim,ydim) hermitian

                    maximize trace(R_win_sum'*rho)
                    subject to
                        trace(rho) == 1;
                        rho >= 0;        
               cvx_end
    
                cval = 1/num_inputs * cvx_optval;            

                if cval > max_cval
                    max_cval = cval;            
                end
            end
        end
    end
end

if num_inputs == 4
    
    for i = 1:num_outputs
        for j = 1:num_outputs
            for k = 1:num_outputs
                for l = 1:num_outputs
                
                    R_win_sum = R{1}{i} + R{2}{j} + R{3}{k} + R{4}{l};

                    cvx_begin sdp quiet 
                        %cvx_precision best; 
                        %#ok<*VUNUS>    % suppress MATLAB warnings for equality checks in CVX
                        %#ok<*EQEFF>    % suppress MATLAB warnings for inequality checks in CVX            
                        variable rho(xdim,ydim) hermitian

                        maximize trace(R_win_sum'*rho)
                        subject to
                            trace(rho) == 1;
                            rho >= 0;        
                   cvx_end

                    cval = 1/num_inputs * cvx_optval;            

                    if cval > max_cval
                        max_cval = cval;            
                    end
                    
                end
            end
        end
    end
end

if num_inputs == 5
    
    for i = 1:num_outputs
        for j = 1:num_outputs
            for k = 1:num_outputs
                for l = 1:num_outputs
                    for m = 1:num_outputs
                
                        R_win_sum = R{1}{i} + R{2}{j} + R{3}{k} + R{4}{l} + R{5}{m};

                        cvx_begin sdp quiet 
                            %cvx_precision best; 
                            %#ok<*VUNUS>    % suppress MATLAB warnings for equality checks in CVX
                            %#ok<*EQEFF>    % suppress MATLAB warnings for inequality checks in CVX            
                            variable rho(xdim,ydim) hermitian

                            maximize trace(R_win_sum'*rho)
                            subject to
                                trace(rho) == 1;
                                rho >= 0;        
                       cvx_end

                        cval = 1/num_inputs * cvx_optval;            

                        if cval > max_cval
                            max_cval = cval;            
                        end
                    end
                end
            end
        end
    end
end

end