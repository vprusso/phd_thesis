%--------------------------------------------------------------------------
% 	NPAR_Search_Monogamy_LB: 	Computes the lower bound of the quantum 
%								value of a refereed non-local game. 
%	This function has 3 required arguments:
%		R: a cell array consisting of the bases for the referee.	
%		REPS: the number of parallel repetitions performed.
%		LVL: The level of the NPAR hierarchy. 
%
% 	BEST = NPAR_SEARCH_MONOGAMY_LB(R,REPS,LVL) is the lower bound quantum 
%	value in a nonlocal game where the referee is allowed to perform a 
%	measurement and is sharing the entangled state with Alice and Bob. This
%	value is obtained via an alternating projection between two SDPs. In 
%   the first SDP, we fix Bob's measurements and optimize over Alice's. In 
%   the second SDP, we fix Alice's measurements and optimize over Bob's.
%
%	This function has two option input arguments:
%		I_MAX:   Times to run the outer loop
%		J_MAX:   Times to run the alternating projection algorithm 
%
%	BEST = NPAR_SEARCH_MONOGAMY_LB(R,REPS,LVL,I_MAX,J_MAX,ALT_MAX) is the 
%   lower bound quantum value as above, but now the precision may be 
%   adjusted by how many times the alternating projection algorithm is run. 
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
function [best, tau, rho, opt_strat_A, opt_strat_B] = MonogamyGameValueLB(R, reps, lvl, varargin)
%#ok<*VUNUS>    % suppress MATLAB warnings for equality checks in CVX
%#ok<*EQEFF>    % suppress MATLAB warnings for inequality checks in CVX 

q0 = [1;0;0]; q1 = [0;1;0]; q2 = [0;0;1]; 

u0 = kron(q1,q2) - kron(q2,q1);
u1 = kron(q0,q2) - kron(q2,q0);
u2 = kron(q0,q1) - kron(q1,q0);

% set optional argument defaults: 
%	I_MAX:   Times to run the outer loop
%	J_MAX:   Times to run the alternating projection algorithm 
[i_max,j_max] = opt_args({ 0, 4 },varargin{:});

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

I = eye(xdim,ydim);

% Setup linear function for monogamy game.
K = zeros(xdim,ydim,num_inputs,num_inputs,num_outputs,num_outputs);
for i = 1:num_inputs
    for j = 1:num_outputs
        K(:,:,i,i,j,j) = R{i}{j};
    end
end
R = K; 
%tau = 1/3*(Normalize(u0*u0') + Normalize(u1*u1') + Normalize(u2*u2'));
%tau = Normalize(1/3*( u0*u0' + u1*u1' + u2*u2' ));
%tau = 1/3 * ( kron(q0,q0)*kron(q0,q0)' + kron(q1,q1)*kron(q1,q1)' + kron(q2,q2)*kron(q2,q2)' );

%cvx_precision best
cvx_solver sedumi
for t = 0:i_max
    best = 0;
        
    for k = 1:j_max
        k
        
        % Generate random bases from the orthogonal colums of randomly 
        % generated unitary matrices.     
        B = zeros(xdim,ydim,num_inputs,num_outputs);
        for y = 1:num_inputs
            U = RandomUnitary(num_outputs);
            for b = 1:num_outputs
                B(:,:,y,b) = U(:,b)*U(:,b)';
            end
        end  

        % Run the actual alternating projection algorithm between
        % the two SDPs. 
        it_diff = 1;
        prev_win = -1;
        while it_diff > 10^-6
            % Optimize over Alice's measurement operators while
            % fixing Bob's. If this is the first iteration, then the 
            % previously randomly generated operators in the outer loop are
            % Bob's. Otherwise, Bob's operators come from running the next
            % SDP.
            cvx_begin sdp quiet                         
                variable rho(xdim^2,ydim^2,num_inputs,num_outputs) hermitian
                variable tau(xdim^2,ydim^2) hermitian

                win = 0;
                for x = 1:num_inputs
                    for y = 1:num_inputs
                        for a = 1:num_outputs
                            for b = 1:num_outputs                               
                                win = win + trace( (kron(R(:,:,x,y,a,b), B(:,:,y,b)))' * rho(:,:,x,a) );                                
                            end
                        end
                    end
                end
                
                
                maximize real(win)

                subject to 
                    
                    % Sum over "a" for all "x". 
                    rho_a_sum = sum(rho,4);
                    for x = 1:num_inputs
                        rho_a_sum(:,:,x) == tau;
                    end
                                                                               
                    % Enforce that tau is a density operator.
                    trace(tau) == 1;
                    tau >= 0;
                                       
                    rho >= 0; 

            cvx_end            
            win = real(1/(num_inputs)*win);
                        
            % Now, optimize over Bob's measurement operators and fix 
            % Alice's operators as those coming from the previous SDP.
            cvx_begin sdp quiet
                variable B(xdim,ydim,num_inputs,num_outputs) hermitian 
                                                             
                win = 0;
                for x = 1:num_inputs
                    for y = 1:num_inputs
                        for a = 1:num_outputs
                            for b = 1:num_outputs                                                               
                                win = win + trace( (kron(R(:,:,x,y,a,b), B(:,:,y,b)))' * rho(:,:,x,a) );
                            end
                        end
                    end
                end     
                
                maximize real(win)
                
                subject to 
                              
                    % Bob's measurements operators must be PSD and sum to I
                    B_b_sum = sum(B,4);
                    for y = 1:num_inputs
                        B_b_sum(:,:,y) == I;
                    end
                    B >= 0;                                                           
                                 
            cvx_end           
            win = real(1/(num_inputs)*win);
            
            it_diff = win - prev_win;
            prev_win = win;
         end
        
        % As the SDPs keep alternating, check if the winning probability
        % becomes any higher. If so, replace with new best.
        if best < win
            
            best = win;
            
            A = zeros(xdim,ydim,num_inputs,num_outputs); 
            for x = 1:num_inputs
                for a = 1:num_outputs
                   A(:,:,x,a) = num_outputs*PartialTrace(rho(:,:,x,a),1); 
                   %A(:,:,x,a) = 1/num_outputs * PartialTrace( tau^(-1/2) * rho(:,:,x,a) * tau^(-1/2) );
                end
            end
            
            opt_strat_A = A;
            opt_strat_B = B;
        end
    
    end;
     
    %best    
end
end