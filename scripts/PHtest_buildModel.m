function res = PHtest_buildModel(prm,simtraj,PHtype)
% res = PHtest_buildModel(prm,simtraj,PHtype)
%
% Generate sythetic FRET state sequences using input parameters and return 
% resulting dwell times and trajectories.
%
% prm: structure with fields:
%  prm.rate: frame rate (per second)
%  prm.N: number of trajectories, set to 0 to use true distributions
%  prm.L: max. trajectory length
%  prm.ndt: minimum number of observed dwell times
%  prm.val: [1-by-J] state values
%  prm.tp: [J-by-J] state transition probability matrix
%  prm.ip: [1-by-J] (opt) initial state probabilities
% simtraj: (1) to return state sequences, (0) to return only dwell times
% PHtype: 1 for discrete PH, 2 for continuous PH
% res: structure with fields:
%  res.dt_gt: {1-by-N} [ndt_gt-by-5] simulated true dwell times (time,state indexes, state values)
%  res.dt_obs: {1-by-N} [ndt_obs-by-3] simulated observed dwell times (time,state indexes, state values)
%  res.seq: (opt) {1-by-N} [L-by-1] simulated state sequences
%
% example:
% >> prm.N = 15;
% >> prm.L = Inf;
% >> prm.ndt = 8;
% >> prm.val = [0.2,0.7,0.7,0.7];
% >> prm.tp = [0,1/300,1/300,1/300; 1/30,0,1/30,1/30; 1/600,1/600,0,1/600; 1/12000,1/12000,1/12000,0];
% >> prm.ip = [0.25,0.25,0.25,0.25];
% >> res = PHtest_buildModel(prm,0)
% 
% res = 
% 
%   struct with fields:
% 
%      dt_gt: {1×15 cell}
%     dt_obs: {1×15 cell}

% CONSTANTS
MAX_TRAJ_LENGTH = 4000;
MIN_PMF = 1E-7;

% initialize results
res = [];

% get number of states and sample size
J = numel(prm.val);

% type of simulation
ndt_based_sim = isinf(prm.L); % based on nb. of observed dwelltimes/traj length
is_mcmc = prm.N>0; % MCMC simulation/use true PDF counts

% transition probabilities
k = prm.tp;
k(~~eye(J)) = 0;
w = k./repmat(sum(k,2),[1,J]);
w(isnan(w)) = 0;

% state lifetimes
tau = 1./sum(k,2);
is_kin_trap = isinf(tau);

% initial state probabilities
isip = isfield(prm,'ip');
if isip % user-defined input
    ip = prm.ip;
else % calculated from state lifetimes
    ip = tau/sum(tau);
end

% observed state aggregates
mu_agg = unique(prm.val);
V = numel(mu_agg);

if ~is_mcmc
    % use true PDF as dwell time count
    ndt = round(prm.ndt/V);
    [T,t,ip] = getPHprm(k,prm.val,PHtype);
    res.dt_obs = cell(1,1);
    res.dt_gt = [];
    for v = 1:V
        [x,P] = calcsimprob(T{v},t{v},ip{v},MIN_PMF,PHtype);
        dt = randsample(x,ndt,true,P);
        P = histcounts(dt,[x,x(end)+1]);
        x = x(P>0);
        P = P(P>0);
        ndat = numel(x);
        res.dt_obs{1} = cat(1,res.dt_obs{1},...
            [x',repmat(v,ndat,1),zeros(ndat,1),P']);
    end
    res.dt_obs{1}(:,1) = res.dt_obs{1}(:,1)/prm.rate;
    return

else
    % initializes results for MCMC simulation
    dt_gt = cell(1,prm.N);
    dt_obs = cell(1,prm.N);
    if simtraj
        seq = cell(1,prm.N);
    end
end

% state aggregate assignment
agg = zeros(1,J);
for j = 1:J
    agg(j) = find(mu_agg==prm.val(j));
end

% pre-determine whether exiting aggregate is allowed
wont_exit_agg = false(1,J);
for j = 1:J
    step = 1;
    j_from = j;
    while step<J
        [~,j_to] = find(w(j_from,:)>0);
        if any(prm.val(j_to)~=prm.val(j))
            break
        end
        j_from = j_to; % prepare next iteration
        step = step+1; % increment exploration step
    end
    if step==(J+1)
        wont_exit_agg(j) = true;
    end
end

% pre-calculate cumulated state probabilities for fast random draws
w_cumsum = cumsum(w,2);
ip_cumsum = cumsum(ip);

% generate dwell times
isTrans = any(any(w>0));
if J>1 && isTrans
    n = 1;
    while n <= prm.N
        % pick a first state
        % state1 = randsample(1:J,1,true,ip); 
        state1 = find(ip_cumsum>=rand,1,'first'); % faster

        l = 0;
        ndt = 0;
        % stes = zeros(1,J); % state mixing coeff in current time bin
        while l<prm.L && ndt<prm.ndt
            if ndt_based_sim && wont_exit_agg(state1)
                % generate state sequence within aggregate only
                while l<MAX_TRAJ_LENGTH
                    % pick a next state
                    % state2 = randsample(1:J,1,true,w(state1,:)); 
                    state2 = find(w_cumsum(state1,:)>=rand,1,'first'); % faster

                    % draw a dwell time
                    % dl = random('exp',tau(state1)); % decimal
                    dl = ceil(-tau(state1)*log(rand())); % integer

                    % append trajectory
                    if simtraj
                        seq{n} = cat(1,seq{n},state1*ones(dl,1));
                    end

                    % append dwell time table
                    [dt_gt{n},dt_obs{n},~] = incrDtTables(dl,state1,state2,...
                        prm.val,ndt,dt_gt{n},dt_obs{n});

                    % prepare next iteration
                    state1 = state2;
                    l = l+dl;
                end
                break
            end
            
            % determine current dwell time
            if is_kin_trap(state1) % kinetic trap
                state2 = state1; % no next state
                if ndt_based_sim
                    % fill trajectory up to max. length
                    dl = MAX_TRAJ_LENGTH-l;
                    if simtraj
                        seq{n} = cat(1,seq{n},state1*ones(dl,1));
                    end

                    % append dwell time table
                    [dt_gt{n},dt_obs{n},~] = incrDtTables(dl,state1,state2,...
                        prm.val,ndt,dt_gt{n},dt_obs{n});
                    break
                else
                    % fill trajectory up to defined length
                    dl = prm.L-l;
                end
            else
                % draw a dwell time
                % dl = random('exp',tau(state1)); % decimal 
                dl = ceil(-tau(state1)*log(rand())); % integer

                % pick a next state
                % state2 = randsample(1:J, 1, true, w(state1,:)); 
                state2 = find(w_cumsum(state1,:)>=rand,1,'first'); % faster
            end
            
            % prevent overflow when trajectory length is finite
            if (l+dl)>prm.L
                dl = prm.L-l;
            end
            
            % append dwell time table
            [dt_gt{n},dt_obs{n},ndt] = incrDtTables(dl,state1,state2,...
                prm.val,ndt,dt_gt{n},dt_obs{n});

            % append trajectory (integer dwell times: no time-averaging of
            % states)
            if simtraj
                seq{n} = cat(1,seq{n},state1*ones(dl,1));
            end
            
            % prepare next iteration
            l = l + dl;
            state1 = state2;

            % % decimal dwell times: handle time-averaging of states
            % if l>0 && sum(stes)<=1
            % 
            %     % the cumulation of the dt generated overflows or 
            %     % reaches the integration time limit
            %     if dl>=(1-sum(stes))
            %         dl = dl - (1 - sum(stes));
            %         l = l + (1 - sum(stes));
            %         stes(state1) = stes(state1) + 1 - sum(stes);
            %         [~,ste] = max(stes);
            %         if simtraj
            %             seq{n} = cat(1,seq{n},ste);
            %         end
            %         stes = zeros(J,1);
            % 
            %     % the cumulation of the dt generated does not reach the
            %     % integration time limit
            %     else
            %         stes(state1) = stes(state1) + dl;
            %         l = l + dl;
            %         state1 = state2;
            %         continue
            %     end
            % end
            % 
            % % append state trajectory with full time bins
            % addl = fix(dl);
            % if addl>=1
            %     if simtraj
            %         seq{n} = cat(1,seq{n},repmat(state1,addl,1));
            %     end
            %     stes = zeros(1,J);
            % end
            % 
            % % keep rest for next iterations
            % fract_end = dl-addl;
            % if fract_end>0
            %     stes(state1) = fract_end;
            % end
            % 
            % % prepare next iteration
            % l = l + fract_end + addl;
            % state1 = state2;
        end

        if simtraj
            for j = 1:J
                seq{n}(seq{n}==j) = agg(j);
            end
        end
        n = n+1;
    end
    
else
    for n = 1:prm.N
        if J > 1
            if any(sum(k,1)==0) || any(sum(k,2)==0)
                % MH 19.12.2019: identify zero sums in rate matrix
                disp(['Simulation aborted: when no transition ',...
                    'is defined (null rates), initial state probabilities',...
                    ' must be pre-defined.']);
                return
            end
            
            % pick a "first state" randomly
            state1 = randsample(1:J, 1, true, sum(ip,1));
        else
            state1 = 1;
        end
        
        if ndt_based_sim
            Ln = MAX_TRAJ_LENGTH;
        else
            Ln = prm.L;
        end

        dt_gt{n} = [Ln state1 NaN prm.val(state1) NaN];
        dt_obs{n} = ...
            [Ln find(prm.val==prm.val(state1)) NaN prm.val(state1) NaN];
        if simtraj
            seq{n} = repmat(state1,Ln,1);
        end
    end
end

% remove dt=0
for n = 1:prm.N
    dt_gt{n} = bindttable(dt_gt{n});
    dt_gt{n}(:,1) = dt_gt{n}(:,1)/prm.rate;
    dt_obs{n} = bindttable(dt_obs{n});
    dt_obs{n}(:,1) = dt_obs{n}(:,1)/prm.rate;
end

% save results
if simtraj
    res.seq = seq;
end
res.dt_gt = dt_gt;
res.dt_obs = dt_obs;


function [dt,prob] = calcsimprob(T,mu,ip,minprob,PHtype)
% dt = [];
% prob = [];
% n = 0;
% while isempty(prob) || prob(end)>minprob
%     n = n+1;
%     dt = cat(2,dt,n);
%     prob = cat(2,prob,ip*(T^(n-1))*mu);
% end

fact = 1.5;
t = 0;
n = 0;
dt = [];
prob = [];
decr = false;
while ~isinf(t) && ~(~isempty(prob) && decr && prob(end)<minprob)
%     if t==0
%         t = 1;
%     else
%         t = t*fact;
%     end

    t = t+1;
    
%     % rounding times prevents zero discrete phase-type PDFs
%     if PHtype==1
%         t = round(t);
%         if ~isempty(dt) && t==dt(end)
%             continue
%         end
%     end
    
    dt = cat(2,dt,t);
    n = n+1;
    switch PHtype
        case 1 % DPH
            if n>1
                CDF1 = calcDPHCDF(T,ip,dt(n-1),eps("double"));
                CDF2 = calcDPHCDF(T,ip,dt(n),eps("double"));
            else
                CDF1 = calcDPHCDF(T,ip,0,eps("double"));
                CDF2 = calcDPHCDF(T,ip,dt(n),eps("double"));
            end
        case 2 % CPH
            if n>1
                CDF1 = calcCPHCDF(T,ip,dt(n-1),eps("double"));
                CDF2 = calcCPHCDF(T,ip,dt(n),eps("double"));
            else
                CDF1 = calcCPHCDF(T,ip,0,eps("double"));
                CDF2 = calcCPHCDF(T,ip,dt(n),eps("double"));
            end
    end
    prob = cat(2,prob,CDF2-CDF1);

    if t>1 && (prob(n)-prob(n-1))<0
        decr = true;
    end
end


function [T,t,ip] = getPHprm(tp,states,PHtype)
val_v = unique(states);
V = numel(val_v);

T = cell(1,V);
t = cell(1,V);
ip = cell(1,V);

tp(~~eye(size(tp))) = 1-sum(tp,2);
for v = 1:V
    idv = states==val_v(v);
    T{v} = tp(idv,idv);
    t{v} = 1-sum(T{v},2);
    ip{v} = estimateiniprob(tp,idv);
    
    if PHtype==2 % CPH
        szTv = size(T{v});
        T{v}(~~eye(szTv)) = T{v}(~~eye(szTv))-1;
    end
end


function [dt_gt,dt_obs,ndt] = incrDtTables(dl,j1,j2,valj,ndt,dt_gt,dt_obs)

valv = unique(valj);

v1 = find(valv==valj(j1));
v2 = find(valv==valj(j2));

dt_gt = cat(1,dt_gt,[dl j1 j2 valj(j1) valj(j2)]);
if ~isempty(dt_obs) && v1==dt_obs(end,2)
    dt_obs(end,1) = dt_obs(end,1)+dl;
    dt_obs(end,3) = v2;
    dt_obs(end,5) = valj(j2);
    if v1~=v2
        ndt = ndt+1;
    end
else
    dt_obs = cat(1,dt_obs,[dl v1 v2 valj(j1) valj(j2)]);
    if v1~=v2
        ndt = ndt+1;
    end
end


function dt = bindttable(dt)
% dt = bindttable(dt)
%
% Round times in dwell times table, remove 0 dwell times and adjust table 
% accordingly.
%
% dt: [ndt-by-3 or 5] dwell time table

dt(dt(:,1)==0,:) = [];
n = 1;
ndt = size(dt,1);
excl = false(1,ndt);
while n<ndt
    m = 1;
    while (n+m)<=size(dt,1) && dt(n+m,2)==dt(n,2)
        dt(n,1) = dt(n,1)+dt(n+m,1);
        excl(n+m) = true;
        m = m+1;
    end
    n = n+m;
end
dt(excl,:) = [];

dt(:,3) = [dt(2:end,2);NaN];
if size(dt,2)==5
    dt(:,5) = [dt(2:end,4);NaN];
end

