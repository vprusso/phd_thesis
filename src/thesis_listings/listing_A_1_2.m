e0 = [1;0]; e1 = [0;1]; ep = [1;1]/sqrt(2); em = [1;-1]/sqrt(2);
psi0_dm = e0*e0'; psi0_dmc = e1*e1';
psi1_dm = ep*ep'; psi1_dmc = em*em';

R00 = psi0_dm/2;
R01 = psi0_dmc/2; 
R10 = psi1_dm/2;
R11 = psi1_dmc/2;  

dim = 9;
A00_B00 = zeros(dim); A01_B01 = zeros(dim); 
A10_B10 = zeros(dim); A11_B11 = zeros(dim); 

% These are the relative positions of these entries as
% indexed by strings in the matrix. 
A00_B00(2,6) = 1; A00_B00(6,2) = 1;
A01_B01(3,7) = 1; A01_B01(7,3) = 1;
A10_B10(4,8) = 1; A10_B10(8,4) = 1;
A11_B11(5,9) = 1; A11_B11(9,5) = 1;

A00_B00 = zeros(dim); A01_B01 = zeros(dim); 

A00_B10 = zeros(dim); A01_B11 = zeros(dim);

A10_B00 = zeros(dim); A11_B01 = zeros(dim); 

A10_B11 = zeros(dim); A11_B10 = zeros(dim); 

A00_B00(2,6) = 1; A00_B00(6,2) = 1;
A00_B10(2,8) = 1; A00_B10(8,2) = 1;

A01_B01(3,7) = 1; A01_B01(7,3) = 1;
A01_B11(3,9) = 1; A01_B11(9,3) = 1;

A10_B00(4,6) = 1; A10_B00(6,4) = 1;
A10_B11(4,9) = 1; A10_B11(9,4) = 1;

A11_B01(5,7) = 1; A11_B01(7,5) = 1;
A11_B10(5,8) = 1; A11_B10(8,5) = 1;

% CHSH ENLG
A = 1/4*(kron(R00, A00_B00) + kron(R01, A01_B01)) + ...
    1/4*(kron(R00, A00_B10) + kron(R01, A01_B11)) + ...
    1/4*(kron(R00, A10_B00) + kron(R01, A11_B01)) + ...
    1/4*(kron(R10, A10_B11) + kron(R11, A11_B10));

cvx_begin sdp
    cvx_precision best
    %#ok<*VUNUS>    % suppress MATLAB warnings for equality checks in CVX
    %#ok<*EQEFF>    % suppress MATLAB warnings for inequality checks in CVX 
    
    % Admissible matrix
    variable M(2*dim,2*dim) hermitian
    
    % Sub-block matrices found in the admissible matrix
    variable M11(dim,dim)
    variable M12(dim,dim)
    
    variable M21(dim,dim)
    variable M22(dim,dim)
            
    M == [ M11 M12;
           M21 M22 ];
       
    maximize trace( A*M )      
    
    subject to 
 
    % Normalization condition:
    M11(1,1) + M22(1,1) == 1;
        
    for i = 1:dim
        for j = 1:dim
            % Ensure commutation relation holds
            %(i.e. [A,B] = 0)            
            M11(i,j) == M11(j,i);
            M12(i,j) == M12(j,i);
            M21(i,j) == M21(j,i);
            M22(i,j) == M22(j,i);
            
            % Enforce operators as projective measurements
            % (i.e. the square of the same operator is found in the top
            % column / row of the diagonal entry). 
            M11(i,i) == M11(1,i);
            M11(i,i) == M11(i,1);
                        
            M12(i,i) == M12(1,i);
            M12(i,i) == M12(i,1);

            M21(i,i) == M21(1,i);
            M21(i,i) == M21(i,1);

            M22(i,i) == M22(1,i);
            M22(i,i) == M22(i,1);
        end
    end  

    % Enforce that projective measurements sum to 1:
    for i = 1:dim
        for j = 1:dim
            if mod(i,2) == 0                
            M11(i,j) + M11(i+1,j) == M11(1,j);
            M12(i,j) + M12(i+1,j) == M12(1,j);
            M21(i,j) + M21(i+1,j) == M21(1,j);
            M22(i,j) + M22(i+1,j) == M22(1,j);
            end
            if mod(j,2) == 0
            M11(i,j) + M11(i,j+1) == M11(i,1);
            M12(i,j) + M12(i,j+1) == M12(i,1);
            M21(i,j) + M21(i,j+1) == M21(i,1);
            M22(i,j) + M22(i,j+1) == M22(i,1);                      
            end
        end
    end

    % Ensure that the matrix is PSD.
    M >= 0;
cvx_end