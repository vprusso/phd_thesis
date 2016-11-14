e0 = [1;0];         e1 = [0;1];
ep = [1;1]/sqrt(2); em = [1;-1]/sqrt(2);

psi0_dm = e0*e0'; psi0_dmc = e1*e1';
psi1_dm = ep*ep'; psi1_dmc = em*em';

lvl = 1;
reps = 1;    
j_max = 4;

xdim = 2;
ydim = 2;

num_inputs = 2;
num_outputs = 2;

I = eye(xdim,ydim);

R = zeros(2,2,2,2,2,2);
R(:,:,1,1,1,1) = psi0_dm/2;
R(:,:,1,1,2,2) = psi0_dmc/2;
R(:,:,2,2,1,1) = psi1_dm/2;
R(:,:,2,2,2,2) = psi1_dmc/2; 

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
            variable rho(xdim^(2*reps),ydim^(2*reps),...
            num_inputs,num_outputs) hermitian
            variable tau(xdim^(2*reps),ydim^(2*reps)) hermitian

            win = 0;
            for x = 1:num_inputs
                for y = 1:num_inputs
                    for a = 1:num_outputs
                        for b = 1:num_outputs                               
                            win = win + ...                            
                            trace( (kron(R(:,:,x,y,a,b), ...
                            B(:,:,y,b)))' * rho(:,:,x,a) );                                
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
        win = real(win);
                   
        % Now, optimize over Bob's measurement operators and fix 
        % Alice's operators as those coming from the previous SDP.
        cvx_begin sdp quiet
                        
            variable B(xdim,ydim,num_inputs,num_outputs) hermitian 
                                                         
            win = 0;
            for x = 1:num_inputs
                for y = 1:num_inputs
                    for a = 1:num_outputs
                        for b = 1:num_outputs                                                               
                            win = win + ...
                             trace( (kron(R(:,:,x,y,a,b), ...
                             B(:,:,y,b)))' * rho(:,:,x,a) );
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
        win = real(win);
        
        it_diff = win - prev_win;
        prev_win = win;
     end
    
    % As the SDPs keep alternating, check if the winning probability
    % becomes any higher. If so, replace with new best.
    if best < win
        
        best = win;

        % take purification of tau
        pur = PartialTrace(tau,2);            

        A = zeros(xdim,ydim,num_inputs,num_outputs); 
        for x = 1:num_inputs
            for a = 1:num_outputs
               A(:,:,x,a) = pur^(-1/2) * PartialTrace(rho(:,:,x,a),2) * pur^(-1/2); 
            end
        end
         
        opt_strat_A = A;
        opt_strat_B = B;
    end

end;
 
best