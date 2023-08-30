function [dphres,expres,fulldphres] = MLDPH_analysis(dest,dt,simprm,prm,...
    fitfull)
% [dphres,expres,fulldphres] = MLDPH_analysis(dest,dt,simprm,prm,fitfull)
%
% dest: export destination folder
% dt: [ndt-by-3] observed dwell time table
% simprm: simulation parameters
% prm: ml-dph parameters
% fitfull: (1) to fit hyperexponential, Erlang and full generator matrices, 
%  (0) for hyperexponential and Elang only

% initialize output
dphres = [];
expres = [];
fulldphres = [];

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

[D,mdlopt,mdl,dthist] = script_findBestModel(dt(:,[1,4,2,3]),Dmax,states,...
    expT,dt_bin,T,false,dest);

% reshape inferrence results
BICres = [];
BICres_sumexp = [];
V = numel(states);
for v = 1:V
    S = size(mdl{v},1);
    for s = 1:S
        Ds = numel(mdl{v}(s,1).pi_fit);
        nfp = sum(sum(mdl{v}(s,1).schm))-1;
        BICs = mdl{v}(s,1).BIC;
        cvg = mdl{v}(s,1).cvg;
        BICres = cat(1,BICres,[v,Ds,s,nfp,BICs,cvg]);
        
        % hyperexponential
        if s==1 || mod(s,2)==0
            BICres_sumexp = cat(1,BICres_sumexp,[v,Ds,s,nfp,BICs,cvg]);
        end
    end
end

% gather results for ground analysis
degen = [];
minBIC = [];
for v = 1:V
    degen = cat(2,degen,repmat(v,[1,D(v)]));
    minBIC = cat(1,minBIC,[mdlopt.BIC(v),D(v)]);
end
states_dph = states(degen);
dphres = {minBIC,BICres,mdlopt,states_dph};

% collect hyperexpoential and full-generator analysis results
mdlopt_sumexp = initmdlstructure(V);
if fitfull
    mdlopt_full = initmdlstructure(V);
end
D_sumexp = zeros(1,V);
for v = 1:V
    
    % determine optimum state degeneracy for hyperexponential fits only
    isvandcvg = BICres_sumexp(:,1)==v & BICres_sumexp(:,6)==1;
    sopt = find(isvandcvg & BICres_sumexp(:,5)==...
        min(BICres_sumexp(isvandcvg,5)));
    sopt = BICres_sumexp(sopt(1),3);
    D_sumexp(v) = size(mdl{v}(sopt,1).schm,1)-2;
    
    % infer hyperexponential parameters on original data
    ffile = [dest,sprintf('_exp_state%iD%i_dphplot',v,D(v))];
    mdl_v = script_inferPH(dthist{v},T,mdl{v}(sopt,1).schm,ffile);
    mdlopt_sumexp.pi_fit{v} = mdl_v.pi_fit;
    mdlopt_sumexp.tp_fit{v} = mdl_v.tp_fit;
    mdlopt_sumexp.schm{v} = mdl_v.schm;
    mdlopt_sumexp.logL(v) = mdl_v.logL;
    mdlopt_sumexp.N(v) = mdl_v.N;
    nfp = sum(sum(mdl_v.schm))-1;
    mdlopt_sumexp.BIC(v) = nfp*log(mdl_v.N)-2*mdl_v.logL;
    
    % infer full-generator parameters on original data
    if fitfull
        ffile = [dest,sprintf('_full_state%iD%i_dphplot',v,D(v))];
        schm_full = ones(D(v)+2);
        schm_full(1,[1,end]) = 0;
        schm_full(:,1) = 0;
        schm_full(end,:) = 0;
        mdl_v = script_inferPH(dthist{v},T,schm_full,ffile);
        mdlopt_full.pi_fit{v} = mdl_v.pi_fit;
        mdlopt_full.tp_fit{v} = mdl_v.tp_fit;
        mdlopt_full.schm{v} = mdl_v.schm;
        mdlopt_full.logL(v) = mdl_v.logL;
        mdlopt_full.N(v) = mdl_v.N;
        nfp = sum(sum(mdl_v.schm))-1;
        mdlopt_full.BIC(v) = nfp*log(mdl_v.N)-2*mdl_v.logL;
    end
    
end

% save computation time
mdlopt_sumexp.t_dphtest = mdlopt.t_dphtest;
if fitfull
    mdlopt_full.t_dphtest = mdlopt.t_dphtest;
end

% gather results for hyperexponential analysis
degen = [];
minBIC = [];
for v = 1:V
    degen = cat(2,degen,repmat(v,[1,D_sumexp(v)]));
    minBIC = cat(1,minBIC,[mdlopt_sumexp.BIC(v),D_sumexp(v)]);
end
states_sumexp = states(degen);
expres = {minBIC,BICres_sumexp,mdlopt_sumexp,states_sumexp};

% gather results for full-generator analysis
if fitfull
    degen = [];
    minBIC = [];
    for v = 1:V
        degen = cat(2,degen,repmat(v,[1,D(v)]));
        minBIC = cat(1,minBIC,[mdlopt_full.BIC(v),D(v)]);
    end
    fulldphres = {minBIC,BICres,mdlopt_full,states_dph};
end


function mdl = initmdlstructure(V)
mdl.pi_fit = cell(1,V);
mdl.tp_fit = cell(1,V);
mdl.schm = cell(1,V);
mdl.logL = zeros(1,V);
mdl.N = zeros(1,V);
mdl.BIC = zeros(1,V);

