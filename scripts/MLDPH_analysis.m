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

[D,mdlopt,mdl] = script_findBestModel(dt(:,[1,4,2,3]),Dmax,states,expT,...
    dt_bin,T,sumexp,dest);

BICres = [];
for v = 1:numel(D)
    for Dtest = 1:size(mdl{v}{1},1)
        nfp = sum(sum(mdl{v}{1}(Dtest,1).schm))-1;
        BICs = mdl{v}{1}(Dtest,1).BIC;
        BICres = cat(1,BICres,[v,Dtest,0,nfp,BICs]);
    end
    S = size(mdl{v}{2},1);
    for s = 1:S
        Ds = numel(mdl{v}{2}(s,1).pi_fit);
        nfp = sum(sum(mdl{v}{2}(s,1).schm))-1;
        BICs = numel(mdl{v}{2}(s,1).BIC);
        BICres = cat(1,BICres,[v,Ds,s,nfp,BICs]);
    end
end

degen = [];
minBIC = [];
for v = 1:numel(states)
    degen = cat(2,degen,repmat(v,[1,D(v)]));
    minBIC = cat(1,minBIC,[mdlopt.BIC(v),D(v)]);
end
states = states(degen);
res = {minBIC,BICres,mdlopt,states};

