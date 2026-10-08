function id = sortStates(stateval,tau)
% id = sortStates(stateval,tau)
%
% Sort input degenerate states according to their lifetimes.
% Sorting is first performed on state values, then on lifetimes.
%
% stateval: [1-by-J] state values
% tau: [1-by-J] state lifetimes
% id: [1-by-J] new state order

id = [];
val = sort(unique(stateval));
for v = 1:numel(val)
    idv = find(stateval==val(v));
    [~,idtau] = sort(tau(idv));
    id = cat(2,id,idv(idtau));
end

