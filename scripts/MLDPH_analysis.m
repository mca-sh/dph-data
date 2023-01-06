function res = MLDPH_analysis(dest,dt,simprm,prm,sumexp)
% MLDPH_analysis(dest,dt,simprm,prm,sumexp)
%
% dest: export destination folder
% dt: [ndt-by-3] observed dwell time table
% simprm: simulation parameters
% prm: ml-dph parameters
% sumexp: (1) to fit sum of exponential, (0) for DPH

% initialize output
res = [];

% get project parameters
expT = 1/simprm.rate;
states = unique(simprm.val);
excl = prm.excl;
dt_bin = prm.bin;
Dmax = prm.Dmax;
T = prm.T;

% get state sequences
N = numel(dt);
dt_new = [];
for n = 1:N
    dt_m = dt{n};
    if size(dt_m,1)<=1 % exclude statics because irreversible transitions give illed distributions
        continue
    end
    
    % remove first and last dwell times
    if excl
        dt_m([1,end],:) = [];
        if size(dt_m,1)<=0
            continue
        end
    end
    dt_new = cat(1,dt_new,[dt_m(:,1:3),n*ones(size(dt_m,1),1)]);
end
dt = dt_new;

if isempty(dt)
    disp(['ML-DPH can not proceed: no dwell times are left after ',...
        'exclusion of static trajectories.']);
    return
end

[D,mdl,BIC,minBIC] = ...
    script_findBestModel(dt(:,[1,4,2,3]),Dmax,states,expT,dt_bin,T,sumexp,...
    dest);

BICres = [];
for D0 = 1:numel(BIC)
    S = size(BIC{D0},2);
    BICres = cat(1,BICres,[repmat(D0,[S,1]),(1:S)',BIC{D0}']);
end

degen = [];
for v = 1:numel(states)
    degen = cat(2,degen,repmat(v,[1,D(v)]));
end
states = states(degen);
res = {minBIC,BICres,mdl,states};

