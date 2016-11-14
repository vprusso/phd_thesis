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

% Level of the extended QC hierarchy corresponds to non-signaling 
lvl = 0;

% Calculate the lower and upper bounds on the BB84 game:
rep_1_val = MonogamyGameValue(R,reps,lvl)

% BB84 game for a single repetition.
reps = 2; 

% Calculate the lower and upper bounds on the BB84 game:
rep_2_val = MonogamyGameValue(R,reps,lvl)