n = 1;
dim = 2^n;

e0 = [1;0];         e1 = [0;1];
ep = [1;1]/sqrt(2); em = [1;-1]/sqrt(2);
eip = (e0 + 1j*e1)/sqrt(2); eim = (e0 - 1j*e1)/sqrt(2);

psi0_dm = e0*e0'; psi0_dmc = e1*e1';
psi1_dm = ep*ep'; psi1_dmc = em*em';
psi2_dm = eip*eip'; psi2_dmc = eim*eim';

P = zeros(2,2,2,2);
%P(:,:,1,1) = (psi0_dm)/2; P(:,:,1,2) = (psi0_dmc)/2;
%P(:,:,2,1) = (psi1_dm)/2; P(:,:,2,2) = (psi1_dmc)/2;

P = zeros(2,2,2,2,2,2);
P(:,:,1,1,1,1) = psi0_dm/2;
P(:,:,1,1,2,2) = psi0_dmc/2;

P(:,:,1,2,1,1) = psi0_dm/2;
P(:,:,1,2,2,2) = psi0_dmc/2;

P(:,:,2,1,1,1) = psi0_dm/2;
P(:,:,2,1,2,2) = psi0_dmc/2;

P(:,:,2,2,1,2) = psi1_dm/2;
P(:,:,2,2,2,1) = psi1_dmc/2;

%%
cvx_begin sdp 
    %#ok<*VUNUS>    % suppress MATLAB warnings for equality checks in CVX
    %#ok<*EQEFF>    % suppress MATLAB warnings for inequality checks in CVX
     
    variable rho(dim,dim,dim,dim,dim,dim) semidefinite
    variable sig(dim,dim,dim,dim) hermitian
    variable xi(dim,dim,dim,dim) hermitian
    variable tau(dim,dim) hermitian 
        
    % construct objective function
    obj_fun = 0;
    for x = 1:dim
        for y = 1:dim
            for a = 1:dim
                for b = 1:dim
                   obj_fun = obj_fun + ip( P(:,:,x,y,a,b), rho(:,:,x,y,a,b) );
               end
           end
       end
    end
    
    maximize  obj_fun
    
    subject to 
    
    rho_b_sum = sum(rho,6);
    for x = 1:dim
        for y = 1:dim
            for a = 1:dim
                rho_b_sum(:,:,x,y,a) == sig(:,:,x,a);
            end
        end
    end
    
    rho_a_sum = sum(rho,5);
    for x = 1:dim
        for y = 1:dim
            for b = 1:dim
                rho_a_sum(:,:,x,y,b) == xi(:,:,y,b);
            end
        end
    end
    
    sig_a_sum = sum(sig,4);
    xi_b_sum = sum(xi,4);
    for x = 1:dim
        sig_a_sum(:,:,x) == tau;
    end
    for y = 1:dim
        xi_b_sum(:,:,y) == tau;
    end
    
    trace(tau) == 1; 
    tau >= 0;          
                        
cvx_end
cvx_optval