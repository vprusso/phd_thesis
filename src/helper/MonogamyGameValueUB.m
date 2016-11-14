%%  MONOGAMYGAMEVALUE    Computes the maximum value of a monogamy game
% 
%   Documentation coming soon (maybe).

function mgval = MonogamyGameValueUB(R,varargin)

    % set optional argument defaults: REPT=1, LVL=1
    [rept,lvl] = opt_args({ 1, 1 },varargin{:});
    
    % Get some basic values and make sure that the input vectors are column
    % vectors.
    num_bases = length(R);
    num_outcomes = length(R{1});
    [xdim,ydim] = size(R{1}{1});
    
    % Now tensor things together if we are doing more than 1 repetition
    if(rept > 1)
        i_ind = zeros(1,rept);
        j_ind = zeros(1,rept);
        
        for i = 1:num_bases^rept
            for j = 1:num_outcomes^rept
                for l = rept:-1:1
                    to_tensor{l} = R{i_ind(l)+1}{j_ind(l)+1};
                end
                newR{i}{j} = Tensor(to_tensor);
                
                j_ind = update_odometer(j_ind,num_outcomes*ones(1,rept));
            end
            i_ind = update_odometer(i_ind,num_bases*ones(1,rept));
        end
        R = newR;
        
        % Recalculate.
        num_bases = length(R);
        num_outcomes = length(R{1});
        [xdim,ydim] = size(R{1}{1});
    end
    
    % Now set up the referee game.
    p = eye(num_bases)/num_bases;
    V = zeros(num_outcomes,num_outcomes,num_bases,num_bases,xdim,ydim);
    for i = 1:num_outcomes
        for j = 1:num_bases
            V(i,i,j,j,:,:) = R{j}{i};
        end
    end

    % Now actually do the hard work; use the modified NPA hierarchy (or
    % whatever other method) to compute the value of this game.
    mgval = RefereeGameValue(p,V,lvl);
end