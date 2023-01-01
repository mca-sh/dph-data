function mdl = script_inferPH(dt,states,expT,dt_bin,n_rs,D,schm,savecurve)
% mdl = script_inferPH(allSchemes,dt,expT,dt_bin,n_rs,J_deg,schm,savecurve)
%
% Trains DPH distributions of specific complexities (in terms of number of degenerated levels) on experimental dwell time histograms (one histogram per state value).
% Returns best fit parameters
%
% dt: [nDt-by-3] dwell times (s), molecule indexes, state values
% states: [1-by-V] state values in dt
% expT: bin time (s)
% dt_bin: binning factor for dwell times prior building histogram
% n_rs: number of restarts
% D: [1-by-V] number of degenerated levels
% schm: {1-by-V}[D(v)-by-D(v)+1] transition schemes to fit
% savecurve: empty or destination folder to save best fit curves
% mdl: structure containing fit DPH parameters
%   mdl.pi_fit: {1-by-V} starting probabilities
%   mdl.tp_fit: {1-by-V} transition probabilities among degenerated states of a same state value
%   mdl.logL: {1-by-V} log likelihoods for best fits
%   mdl.N: [1-by-V] number of data

% defaults
PH_type = 1;% 1 for discrete, 2 for continuous

% initialize output
mdl = struct;

% collect DPH fit curve export
saveit = ~isempty(savecurve);

% gather degenerated state indexes
V = numel(states);
n = 0;
degen = cell(1,numel(D));
for v = 1:numel(D)
    degen{v} = [];
    for j = 1:D(v)
        n = n+1;
        degen{v} = [degen{v} n];
    end
end

% bin dwell times
expT_bin = dt_bin*expT;
dt(:,1) = round(dt(:,1)/expT_bin);

% get oberved dwell time histograms
P = cell(1,V);
x = P;
for v = 1:V
    dt_z = dt(dt(:,3)==v,1);
    if isempty(dt_z)
        continue
    end
    edg = 0.5:(max(dt_z)+0.5);
    x{v} = mean([edg(2:end);edg(1:end-1)],1);
    P{v} = histcounts(dt_z,edg);
end

% calculate phase-type complementary CDF for each state value z
pi_fit = cell(1,V);
tp_fit = pi_fit;
logL = -Inf(1,V);
nDat = zeros(1,V);
if saveit
    P_fit = pi_fit;
    w_fit = pi_fit;
    tau_fit = pi_fit;
end
for v = 1:V
    if isempty(P{v})
        pi_fit{v} = 1;
        tp_fit{v} = [1,0];
        fprintf('state %i/%i: trap state (irreversible transition)\n',v,V);
        continue
    end
    
    v_e = ones(D(v),1);

    pi_fit{v} = NaN(1,D(v));
    tp_fit{v} = NaN(D(v),D(v)+1);
    if saveit
        L = numel(x{v});
        P_fit{v} = zeros(1,L);
        tau_fit{v} = NaN(D(v),1);
        w_fit{v} = NaN(D(v)+1,D(v)+1);
    end
    
    incl = (x{v}>0 & P{v}>0);

    % generate random PH parameters
    a_fit = [];
    T_fit = [];
    for rs = 1:n_rs
        
        fprintf('state %i/%i, restart %i/%i: ',v,V,rs,n_rs);

        % use random starting guess
        a_start = ones(1,D(v));
        a_start = a_start/sum(a_start);

        if PH_type==2 % continuous PH
            w0 = rand(D(v),D(v)+1);
            w0([~~eye(D(v)) false(D(v),1)]) = 0;
            w0 = w0./repmat(sum(w0,2),[1,D(v)+1]);

            r0 = rand(D(v),1);

            T_start = w0.*repmat(r0,[1,D(v)+1]);
            t = T_start(:,end);
            T_start = T_start(:,1:D(v));
            T_start(~~eye(D(v))) = -(sum(T_start,2)+t);
        end

        if PH_type==1 % discrete PH
            tp0 = rand(D(v),D(v)+1);
            tp0(~schm{v}) = 0;
            tp0(~~eye(size(tp0))) = 10;
            tp0 = tp0./repmat(sum(tp0,2),[1,D(v)+1]);

            T_start = tp0(:,1:D(v));
        end

        % train a PH model on experimental CDF
        [a_res,T_res,logL_res] = ...
            trainPH(a_start,T_start,[x{v}(incl);P{v}(incl)]);
        if isempty(a_res) || isempty(T_res)
            continue
        end
        if logL_res>logL(v)
            logL(v) = logL_res;
            a_fit = a_res;
            T_fit = T_res;
        end
    end
    if isempty(a_fit) || isempty(T_fit)
        if saveit
            P_fit{v}(1,:) = 0;
        end
        continue
    end
    
    if saveit
        if PH_type==1 % discrete PH
            t = v_e-T_fit*v_e;
        end
        for l = 1:L
            if PH_type==1 % discrete PH
                P_fit{v}(1,l) = a_fit*(T_fit^(x{v}(l)-1))*t;
            else% continuous PH
                P_fit{v}(1,l) = a_fit*expm(T_fit*x{v}(l))*v_e;
            end
        end
        
        % export hitstogram and DPH fit curve
        dat = [x{v}',P{v}',P_fit{v}'];
        ffile = [savecurve,sprintf('_state%iD%i_dphplot',v,D(v))];
        save(ffile,'dat','-ascii')
    end

    % get parameters from trained PH model
    pi_fit{v} = a_fit;
    if PH_type==2 % continuous PH
        r_v = -diag(T_fit)/dt_bin;
        w = T_fit./repmat(r_v,1,D(v));
        w(~~eye(D(v))) = 0;
        w = cat(2,w,1-sum(w,2));
        tp_fit{v} = w.*repmat(r_v,[1,D(v)+1]);
        tp_fit{v}(~~eye(D(v))) = 1-sum(tp_fit{v},2);
    else % discrete PH
        t = v_e-T_fit*v_e;
        tp_fit{v} = [T_fit,t];
    end
    nDat(v) = sum(P{v});
end
mdl.pi_fit = pi_fit;
mdl.tp_fit = tp_fit;
mdl.logL = logL;
mdl.N = nDat;
mdl.schm = schm;


