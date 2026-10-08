function [T,t,ip,schm] = collectsdegenparam(prm)
% [T,t,ip,schm] = collectsdegenparam(prm)
%
% Collects DPH parameters of each state aggregate from simulation 
% parameters.
%
% prm: structure containing simulation parameters with fields:
%   prm.ip: [1-by-J] initial state probabilities
%   prm.tp: [J-by-J] transition probability matrix
%   prm.val:[1-by-J] FRET state values
% T: {1-by-V}[Dv-by-Dv] transition probabilities between agglomerate states
% t: {1-by-V}[Dv-by-1] exit probability vector
% ip: {1-by-V}[1-by-Dv] initial state probability vector
% schm: {1-by-V}[(Dv+1)-by-(Dv+1)] transition scheme with state1 being the
%       initiating AND absorbing state

val_v = unique(prm.val);
V = numel(val_v);

T = cell(1,V);
t = cell(1,V);
ip = cell(1,V);
schm = cell(1,V);

% propup transition probability matrix
tp0 = prm.tp;
tp0(~~eye(size(tp0))) = 0;
tp0(~~eye(size(tp0))) = 1-sum(tp0,2);

ip0 = prm.ip;

Dv = zeros(1,V);
for v = 1:V
    idv = prm.val==val_v(v);
    Dv(v) = sum(idv);

    T{v} = tp0(idv,idv);
    t{v} = 1-sum(T{v},2);
    
    ip{v} = estimateiniprob(tp0,ip0,idv);
    schm{v} = double([[0,ip{v}];[t{v},tp0(idv,idv)]]>0);
end

