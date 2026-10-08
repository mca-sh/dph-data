function res = PHtest_MLPHanalysis(dest,dt,simprm,prm,verbose)
% res = PHtest_MLPHanalysis(dest,dt,simprm,prm,verbose)
%
% dest: export destination folder
% dt: {1-by-N}[ndt-by-(3 or 4)] observed dwell time (and counts) table
% simprm: simulation parameters with fields:
%   simprm.rate: experiment sampling rate (in s-1)
%   simprm.val: [1-by-J] state values
% prm: ml-ph parameters with fields:
%   prm.excl: 1 to exclude first & last dwell times of each trajectory, 0
%             otherwise.
%   prm.bin: bin size
%   prm.Dmin: min. number of degenerate states to fit
%   prm.Dmax: max. number of degenerate states to fit
%   prm.T: number of model initialization
%   prm.dlmin: convergence criterion on log-likelihood ratio
%   prm.dprmmin: convergence criterion on parameters
%   prm.M_max: maximum nb of EM iterations
%   prm.PHtype: 1 to fit discrete distribution, 2 for continuous
%   prm.applyrules: apply additional selection rules on inferred model 
%                   (distribution validity and state doublons)
%   prm.model_selection: 'cross-validation' or 'BIC' (default)
%   prm.calcmode:-2: EM-EMM with positive coefficients
%                -1: EM-EMM with negative coefficients
%                 0: test canonical connectivities (+10 for searching for 
%                    min. complexity)
%                 1: test all connectivities (+10 for searching for min. 
%                    complexity)
%                 2: test "uncoupled" and "irreversible loop" 
%                    connectivities (+10 for searching for min. complexity)
%                 3: test "uncoupled" connectivity (+10 for searching for 
%                    min. complexity)
%                 4: test "coupled"  (+10 for searching for min. 
%                    complexity)
%                 5: test "uncoupled" and "generalized coxian" connectivity 
%                    (+10 for searching for min. complexity)
%                 6: test "acyclic" connectivity  (+10 for searching for 
%                    min. complexity).
%                 7: test "acyclic" connectivity with different state 
%                    initiations.
%                10: iEMM from Hines et. al. 2015
% verbose: true to show logs, false to mute

% CONSTANTS
NITERS_MAX = 20000; % [iEMM] max. nb. of Gibbs sampling iterations

% initialize output
res = [];

% get project parameters
rate = simprm.rate;
states = unique(simprm.val);
excl = prm.excl;
dt_bin = prm.bin;
Dmax = prm.Dmax;
Dmin = prm.Dmin;
applyrules = prm.applyrules;
calcmode = prm.calcmode;
analysismethod = PHtest_getmethodfromcalcmode(calcmode);
if strcmp(analysismethod,'mlph')
    model_selection = prm.model_selection;
else
    model_selection = 'BIC';
end

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
    dt_new = cat(1,dt_new,dt_m);
end
dt = dt_new;
if isempty(dt) 
    if verbose
        disp(['ML-PH can not proceed: no dwell times are left after ',...
            'exclusion of static trajectories.']);
    end
    return
end

switch analysismethod
    case 'mlph'
        PHtype = prm.PHtype;
        T = prm.T;
        M_max = prm.M_max;
        dlmin = prm.dlmin;
        dprmmin = prm.dprmmin;
        [D,mdlopt,mdl,~] = script_findBestModel(dt,PHtype,Dmin,Dmax,states,...
            rate,dt_bin,T,dlmin,dprmmin,M_max,dest,calcmode,model_selection,...
            applyrules,verbose);

    case 'emexp'
        if calcmode==-1
            negcoeff = true;
        else
            negcoeff = false;
        end
        T = prm.T;
        M_max = prm.M_max;
        dlmin = prm.dlmin;
        dtaumin = prm.dtaumin;

        val = unique(dt(:,2));
        V = numel(val);
        D = zeros(1,V);
        mdl = cell(1,V);
        mdlopt.D = zeros(1,V);
        mdlopt.a = cell(1,V);
        mdlopt.tau = cell(1,V);
        mdlopt.BIC = Inf(1,V);
        mdlopt.logL = -Inf(1,V);
        mdlopt.cvg = false(1,V);
        t = tic;
        for v = 1:V
            dt_v = rate*dt(dt(:,2)==val(v),1)';
            if isequal(M_max,'auto')
                M_max_v = 2*length(dt_v);
            else
                M_max_v = M_max;
            end
            [tau,a,logL,BIC,D(v),tcomp] = BIC_exp(dt_v,Dmin(v),Dmax(v),...
                negcoeff,T,M_max_v,dlmin,dtaumin,verbose);
            for Ds = Dmin(v):Dmax(v)
                mdl_s = [];
                mdl_s.schm = zeros(Ds+2);
                mdl_s.D = Ds;
                mdl_s.a = a{Ds};
                mdl_s.tau = tau{Ds};
                mdl_s.BIC = BIC(Ds);
                mdl_s.logL = logL(Ds);
                mdl_s.cvg = ~isinf(logL(Ds)) & ~isnan(logL(Ds));
                mdl_s.t_emexp = tcomp;
                mdl{v} = cat(1,mdl{v},mdl_s);
                if Ds==D(v)
                    mdlopt.schm{v} = zeros(Ds+1,Ds+2);
                    mdlopt.D(v) = Ds;
                    mdlopt.a{v}= a{Ds};
                    mdlopt.tau{v} = tau{Ds};
                    mdlopt.BIC(v) = BIC(Ds);
                    mdlopt.logL(v) = logL(Ds);
                    mdlopt.cvg(v) = ~isinf(logL(Ds)) & ~isnan(logL(Ds));
                end
            end
        end
        mdlopt.t_emexp = toc(t);

    case 'iemm'
        Kinit = prm.Kinit;
        niters = prm.niters;
        thin = prm.thin;
        hypers = prm.hypers;

        val = unique(dt(:,2));
        V = numel(val);
        D = zeros(1,V);
        mdl = cell(1,V);
        mdlopt.D = zeros(1,V);
        mdlopt.a = cell(1,V);
        mdlopt.tau = cell(1,V);
        mdlopt.BIC = Inf(1,V);
        mdlopt.logL = -Inf(1,V);
        mdlopt.cvg = false(1,V);
        mdlopt.Diter = cell(1,V);
        t = tic;
        for v = 1:V
            t2 = tic;
            dt_v = rate*dt(dt(:,2)==val(v),1)';
            if isequal(niters,'auto')
                niters_v = min([2*length(dt_v),NITERS_MAX]);
            else
                niters_v = niters;
            end
            [Diter,qiter,aiter] = iEMM(dt_v,hypers,Kinit,niters_v,thin);
            tcomp = toc(t2);
            D(v) = mode(Diter);
            for Ds = Dmin(v):Dmax(v)
                iter = find(Diter==Ds);
                mdl_s = [];
                mdl_s.schm = zeros(Ds+2);
                mdl_s.D = Ds;
                mdl_s.t_iemm = tcomp;
                if isempty(iter)
                    mdl_s.a = [];
                    mdl_s.tau = [];
                    mdl_s.BIC = Inf;
                    mdl_s.logL = -Inf;
                    mdl_s.cvg = false;
                else
                    mdl_s.a = mean(cell2mat(aiter(iter)));
                    mdl_s.tau = 1./mean(cell2mat(qiter(iter)));
                    mdl_s.logL = nnz(Diter==Ds)/numel(Diter);
                    mdl_s.BIC = 1/mdl_s.logL;
                    mdl_s.cvg = ~isinf(mdl_s.logL) & ~isnan(mdl_s.logL);
                end
                mdl{v} = cat(1,mdl{v},mdl_s);
                if Ds==D(v)
                    mdlopt.schm{v} = zeros(Ds+1,Ds+2);
                    mdlopt.D(v) = Ds;
                    mdlopt.a{v}= mdl_s.a;
                    mdlopt.tau{v} = mdl_s.tau;
                    mdlopt.BIC(v) = mdl_s.BIC;
                    mdlopt.logL(v) = mdl_s.logL;
                    mdlopt.cvg(v) = mdl_s.cvg;
                end
            end
            mdlopt.Diter{v} = Diter;
        end
        mdlopt.t_iemm = toc(t);
end

% reshape inferrence results
BICres = [];
V = numel(states);
for v = 1:V
    S = size(mdl{v},1);
    for s = 1:S
        switch analysismethod
            case 'mlph'
                Ds = numel(mdl{v}(s,1).pi_fit);
                nfp = sum(sum(mdl{v}(s,1).schm))-1;
                BICs = mdl{v}(s,1).BIC;
                logL_v = mdl{v}(s,1).logL_validation;
                cvg = mdl{v}(s,1).cvg;
                BICres = cat(1,BICres,[v,Ds,s,nfp,BICs,logL_v,cvg]);
            case 'emexp'
                Ds = mdl{v}(s,1).D;
                nfp = 2*Ds-1;
                BICs = mdl{v}(s,1).BIC;
                cvg = mdl{v}(s,1).cvg;
                BICres = cat(1,BICres,[v,Ds,s,nfp,BICs,cvg]);
        end
    end
end

% gather results for ground analysis
degen = [];
minBIC = [];
for v = 1:V
    degen = cat(2,degen,repmat(v,[1,D(v)]));
    switch model_selection
        case 'cross-validation'
            minBIC = cat(1,minBIC,[mdlopt.logL_validation(v),D(v)]);
        otherwise
            minBIC = cat(1,minBIC,[mdlopt.BIC(v),D(v)]);
    end
end
states_ph = states(degen);
res = {minBIC,BICres,mdlopt,states_ph,mdl};

