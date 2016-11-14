% Create the BB84 basis.
e0 = [1;0]; e1 = [0;1];
ep = [1;1]/sqrt(2); em = [1;-1]/sqrt(2); 

psi0 = e0*e0'; psi1 = e1*e1';
psip = ep*ep'; psim = em*em'; 

% Referee's first basis: {|0><0|, |1><1|}
R{1} = {psi0,psi1};

% Referee's second basis: {|+><+|, |-><-|}
R{2} = {psip,psim};

% BB84 game for a single repetition.
reps = 1; 

% Level of the extended QC hierarchy 
lvl = 1;

% Calculate the lower and upper bounds on the BB84 game:
%   cos^2(pi/8) \approx 0.8536
lb = MonogamyGameValueLB(R,reps,lvl)
ub = MonogamyGameValueUB(R,reps,lvl)