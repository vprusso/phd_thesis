% Number of inputs and outputs
nin = 4;
nout = 3;

% Create the mutually unbiased bases consisting of 4-inputs and 3-outputs.
m = MUB(nout); 
R = {};
for i = 1:nin
    for j = 1:nout
        R{i}{j} = m{i}(:,j) * m{i}(:,j)';
    end
end

% Number of repetitions of the game. 
reps = 1; 

% Level of the extended QC hierarchy.
lvl = 1;

% Calculate the lower and upper bounds on the quantum value of
% the mutually unbiased basis game:
lb = MonogamyGameValueLB(R,reps,lvl)
ub = MonogamyGameValueUB(R,reps,lvl)