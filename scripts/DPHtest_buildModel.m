function res = DPHtest_buildModel(prm,simtraj)
% res = DPHtest_buildModel(prm,simtraj)
%
% Generate sythetic FRET state sequences using input parameters and return 
% resulting dwell times and trajectories.
%
% prm: structure with fields:
%  prm.rate: frame rate (per second)
%  prm.N: number of trajectories
%  prm.L: max. trajectory length
%  prm.ndt: minimum number of observed dwell times
%  prm.val: [1-by-J] state values
%  prm.tp: [J-by-J] state transition probability matrix
%  prm.ip: [1-by-J] (opt) initial state probabilities
% simtraj: (1) to return state sequences, (0) to return only dwell times
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
% >> res = DPHtest_buildModel(prm,0)
% 
% res = 
% 
%   struct with fields:
% 
%      dt_gt: {1×15 cell}
%     dt_obs: {1×15 cell}

% defaults
res = [];
Lmax = 4000;

% get number of states and sample size
J = numel(prm.val);

% transition probabilities
prm.tp(~~eye(J)) = 0;
w = prm.tp./repmat(sum(prm.tp,2),[1,J]);
w(isnan(w)) = 0;

% initial state probabilities
isip = isfield(prm,'ip');
tau = 1./sum(prm.tp,2);
if isip % user-defined input
    ip = prm.ip;
else % calculated from state lifetimes
    ip = tau/sum(tau);
end

% initializes results
if simtraj
    if ~isinf(prm.L)
        discr_seq = -ones(prm.L,prm.N);
    else
        discr_seq = [];
    end
end
dt_gt = cell(1,prm.N);
dt_obs = cell(1,prm.N);
if simtraj
    seq = cell(1,prm.N);
end

% get state value indexes
valv = unique(prm.val);
idv = zeros(1,J);
for j = 1:J
    idv(j) = find(valv==prm.val(j));
end

% generate dwell times
isTrans = any(any(w>0));
if J>1 && isTrans
    n = 1;
    while n <= prm.N
        if n <= prm.N
            stes = zeros(1,J);
            l = 0;
            ndt = 0;
            
            state1 = randsample(1:J, 1, true, ip); % pick a first state
            
            while l<prm.L && ndt<prm.ndt
                
                % generate sequence of degenerate states only
                if isinf(prm.L) && isequal(unique(prm.val(w(state1,:)>0)),...
                        prm.val(state1))
                    while l<Lmax
                        state2 = randsample(1:J, 1, true, w(state1,:)); 
                        dl = round(random('exp',tau(state1)));
                        if simtraj
                            seq{n} = cat(1,seq{n},repmat(state1,dl,1));
                        end
                        [dt_gt{n},dt_obs{n},~] = incrDtTables(dl,state1,...
                            state2,prm.val,ndt,dt_gt{n},dt_obs{n});
                        state1 = state2;
                        l = l+dl;
                    end
                    break
                end
                
                % pick a next state
                if sum(w(state1,:))==0
                    state2 = state1;
                else
                    state2 = randsample(1:J, 1, true, w(state1,:)); 
                end
                
                % dwell in number of time bins
                if isinf(tau) % kinetic trap
                    if isinf(prm.L)
                        dl = Lmax-l;
                        if simtraj
                            seq{n} = cat(1,seq{n},repmat(state1,dl,1));
                        end
                        [dt_gt{n},dt_obs{n},~] = incrDtTables(dl,state1,...
                            state2,prm.val,ndt,dt_gt{n},dt_obs{n});
                        break
                    else
                        dl = prm.L-l;
                    end
                else
                    dl = round(random('exp',tau(state1)));
                end
                
                % truncate dwell if larger than remaining time
                if (l+dl)>prm.L
                    dl = prm.L-l;
                end
                
                % update dwell time tables
                [dt_gt{n},dt_obs{n},ndt] = incrDtTables(dl,state1,...
                    state2,prm.val,ndt,dt_gt{n},dt_obs{n});

                if l>0 && sum(stes)<=1
                    
                    % the cumulation of the dt generated overflows or 
                    % reaches the integration time limit
                    if dl>=(1-sum(stes))
                        dl = dl - (1 - sum(stes));
                        l = l + (1 - sum(stes));
                        stes(state1) = stes(state1) + 1 - sum(stes);
                        [~,ste] = max(stes);
                        if simtraj
                            seq{n} = cat(1,seq{n},ste);
                        end
                        stes = zeros(J,1);
                    
                    % the cumulation of the dt generated does not reach the
                    % integration time limit
                    else
                        stes(state1) = stes(state1) + dl;
                        l = l + dl;
                        state1 = state2;
                        continue
                    end
                end

                addl = fix(dl);
                if addl>=1
                    if simtraj
                        seq{n} = cat(1,seq{n},repmat(state1,addl,1));
                    end
                    stes = zeros(1,J);
                end
                
                fract_end = dl-addl;
                % the cumulation of the dt generated overflows the 
                % integration time limit
                if fract_end>0
                    stes(state1) = fract_end;
                end
                
                l = l + fract_end + addl;
                state1 = state2;
            end
            if simtraj
                for j = 1:J
                    seq{n}(seq{n}==j) = idv(j);
                end
            end
            n = n+1;
        end
    end
    
else
    for n = 1:prm.N
        if J > 1
            if any(sum(prm.tp,1)==0) || any(sum(prm.tp,2)==0)
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
        
        if isinf(prm.L)
            Ln = Lmax;
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


function [dt_gt,dt_obs,ndt] = incrDtTables(dl,j1,j2,valj,ndt,dt_gt,dt_obs)

valv = unique(valj);

v1 = find(valv==valj(j1));
v2 = find(valv==valj(j2));

dt_gt = cat(1,dt_gt,[dl j1 j2 valj(j1) valj(j2)]);
if isempty(dt_obs)
    dt_obs = [dl v1 v2 valj(j1) valj(j2)];
    if v1~=v2
        ndt = ndt+1;
    end
elseif valj(j1)==dt_obs(end,2)
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
