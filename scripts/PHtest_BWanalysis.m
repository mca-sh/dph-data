function res = PHtest_BWanalysis(simprm,dphprm,dt,seq,dphres,T,dphguess,...
    verbose)
% res = PHtest_BWanalysis(simprm,dphprm,dt,seq,dphres,T,dphguess,verbose)
%
% Perform Baum-Welch inference on simulated state sequences using ML-PH
% state configuration.
%
% simprm: structure with fields:
%  simprm.rate: simulation frame rate (per second)
%  simprm.val: simulation observed state values
% dphprm: structure with fields:
%  dphprm.excl: (1) to exclude first and last dwell time of each trajectory, (0) otherwise
% dt: {1-by-N} [ndt-by-3] simulated dwell time tables
% seq: {1-by-N} simulated FRET state sequences
% dphres: {1-by-4} results of ML-PH run
% T: number of TPM initializations in BW inferences
% dphguess: (1) to use ML-PH outcome in initial TPM guess, (0) to use only random initial TPM guesses
% verbose: true to show logs, false to mute
% res: {1-by-4} BW results

% defaults
tpmin = 1E-5; % minimum transition probability

% get project parameters
minBIC = dphres{1};
mdl = dphres{3};
states = dphres{4};

% initializes time count
t0 = tic;

% get state sequences
dt_new = [];
N = numel(dt);
Ls = zeros(1,N);
exclmols = false(1,N);
for n = 1:N
    dt_m = dt{n};
    
    % remove first and last dwell times
    if dphprm.excl
        dt_m([1,end],:) = [];
        if size(dt_m,1)<=0
            exclmols(n) = true;
            continue
        end
    end
    dt_new = cat(1,dt_new,[dt_m(:,1:3),n*ones(size(dt_m,1),1)]);
    Ls(n) = round(sum(dt_m(:,1))*simprm.rate);
end
dt = dt_new;
seq(exclmols) = [];

% get relative number of transitions
stateval = unique(simprm.val);
V = numel(stateval);
for v = 1:V
    dt(dt(:,2)==v,2) = stateval(v);
    dt(dt(:,3)==v,3) = stateval(v);
end
clstPop = zeros(V);
for v1 = 1:V
    for v2 = 1:V
        if v1==v2
            continue
        end
        clstPop(v1,v2) = size(dt(dt(:,2)==stateval(v1) & ...
            dt(:,3)==stateval(v2),:),1);
    end
end
clstPop = clstPop/sum(sum(clstPop));

% calculate initial transition prob matrix
D = zeros(1,V);
for v = 1:V
    D(v) = minBIC(v,2);
end
J = sum(D);
if dphguess
    tp0 = zeros(J);
    j1 = 0;
    for v1 = 1:V
        tp0((j1+1):(j1+D(v1)),(j1+1):(j1+D(v1))) = ...
            mdl.tp_fit{v1}(:,1:end-1);
        p_exit = mdl.tp_fit{v1}(:,end);
        j2 = 0;
        for v2 = 1:V
            if v1~=v2
                tp0((j1+1):(j1+D(v1)),(j2+1):(j2+D(v2))) = ...
                    repmat(p_exit,[1,D(v2)]).*...
                    repmat(clstPop(v1,v2),[D(v1),D(v2)])/D(v2);
            end
            j2 = j2+D(v2);
        end
        j1 = j1+D(v1);
    end
else
    tp0 = [];
end

% convert state sequences to uint8 to save memory
for n = 1:numel(seq)
    seq{n} = uint8(seq{n});
end

expPrm.expT = 1/simprm.rate;
expPrm.Ls = Ls;
expPrm.dt = dt(:,[1,4,end-1,end]);
expPrm.excl = dphprm.excl;
expPrm.seq = seq;
schm = true(J);
tp = ones(J);
iter = 1;
eps = 1/(sum(expPrm.Ls)-2*numel(expPrm.Ls));

if verbose
    fprintf('Perform BW iterative search...\n');
end
nb = 0;
nb2 = 0;
while iter==1 || any(tp(schm)<eps)
    if verbose
        nb = dispProgress(sprintf('iter %i:',iter),nb+nb2);
    end
    schm = schm & tp>=tpmin;
    [tp,err,ip,simdat,nb2] = optimizeProbMat(states, expPrm, tp0, T, schm,...
        verbose); % transition prob
    iter = iter+1;
end
if nb>0
    fprintf('\n');
end

res = {tp,err,ip,simdat,toc(t0)};


