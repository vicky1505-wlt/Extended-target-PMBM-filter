function [assignments,costs]= mbestwrap_updt_custom(P0,m,w)

if m==0
    assignments= [];
    costs= [];
    return;
end

n = size(P0,2);

l = length(w);
blk1 = Inf*ones(l);
blk1(logical(eye(l))) = -log(w);

P0 = [P0 blk1];

% Make costs non-negative (required by 'assignmentoptimal')
x = min(min(P0));
P0 = P0 - x;

% Murty
[assignments, costs] = murty_custom(P0,m);

% Restore correct costs to assignments
costs = costs + (x.*sum(assignments>0,2))';

assignments(assignments>n) = 0;

end
