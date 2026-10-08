function [T,t,ip,schm,nb] = calc_absmm_from_hmm(ip0,tp0,states)
% [T,t,ip,schm,nb] = calc_absmm_from_hmm(ip0,tp0,states)
%
% Calculate parameters of absorbing Markov models for each satte aggregate 
% from hidden Markov model parameters.
%
% ip0: [1-by-J] HMM initial state probabilities
% tp0: [J-by-J] HMM state transition probabilities
% states: [1-by-J] state values
%
% T: {1-by-nAgg}[Jagg-by-Jagg] state transition probabilities within the
%    aggregate.
% t: {1-by-nAgg}[Jagg-by-1] state absorption probabilities
% schm: {1-by-nAgg}[(Jagg+1)-by-(Jagg+1)] transition scheme where 0 stands
%       for forbidden and 1 for allowed transitions. State 1 is the
%       initiation AND absorption state.
% nb: number of charcaters printed in prompt

% defaults
nb = 0;

tp0(~~eye(size(tp0))) = 0;
tp0(~~eye(size(tp0))) = 1-sum(tp0,2);

val = unique(states);
V = numel(val);
T = cell(1,V);
t = cell(1,V);
schm = cell(1,V);
ip = cell(1,V);
for v = 1:V
    idv = find(states==val(v));
    ip{v} = estimateiniprob(tp0,ip0,idv);
    T{v} = tp0(idv,idv);
    t{v} = sum(tp0(idv,states~=val(v)),2);
    % Dv = numel(idv);
    % schm{v} = [0,ip{v},0;zeros(Dv,1),T{v},t{v};zeros(1,Dv+2)];
    % if any(isnan(schm{v}(:)))
    %     error('NaN transition probabilities.');
    % end
    % schm{v} = ~~schm{v};
    schm{v} = double([[0,ip{v}];[t{v},tp0(idv,idv)]]>0);
end