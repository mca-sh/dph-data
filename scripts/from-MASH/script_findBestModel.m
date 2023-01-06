function [degen,mdl,BIC,minBIC] = script_findBestModel(dt,Dmax,states,...
    expT,dt_bin,T,sumexp,savecurve)
% [degen,mdl,BIC,minBIC] = script_findBestModel(dt,Dmax,states,expT,dt_bin,T,sumexp,savecurve)
%
% Import dwell times from .clst file
% Find and return most sufficient model complexities (in terms of number of degenerated levels) for each state value
% Plot PH fits and BIC results
%
% dt: [nDt-by-3] dwell times (s), molecule indexes, state values
% Dmax: maximum number of degenerated levels
% states: [1-by-V] state values in dt
% expT: bin time (s)
% dt_bin: binning factor for dwell times prior building histogram
% T: number of restart
% sumexp: (1) to fit a sume of exponential (0) to fit DPH
% savecurve: empty or destination folder to save best fit curves
% degen: [1-by-V] most sufficient model complexity (number of degenerated levels per state value)
% mdl: structures containing best DPH fit parameters for the most sufficient model
%   mdl.pi_fit: {1-by-V} starting probabilities
%   mdl.tp_fit: {1-by-V} transition probabilities among degenerated states of a same state value
%   mdl.logL: {1-by-V} log likelihoods for best fits
%   mdl.N: [1-by-V] number of data
%   mdl.schm: transition scheme
% BIC: {1-by-Dmax}[V-by-S] BIC for all complexities
% minBIC = [V-by-3] min-BIC (1st column), most sufficient state degeneracy
%   (2nd) and transition scheme index (3rd)

% initialize computation time
t_comp = tic;

% Get optimum DPHs for each model complexity
disp('Train DPH distributions on binned dwelltime histograms...')
V = numel(states);
mdl = cell(1,Dmax);
logL = cell(1,Dmax);
for D = 1:Dmax
    fprintf('for %i degenerated states:\n',D);
    for v = 1:V
        if sumexp
            schm0 = [false(D),true(D,1)];
        else
            schm0 = ~eye(D,D+1);
%             schm0 = addExitProb2TransSchemes(getTransSchemes(D));
        end
    end
    S = size(schm0,3);
    logL{D} = -Inf(V,S);
    mdl{D} = struct([]);
    for s = 1:S
        schm = cell(1,V);
        for v = 1:V
            schm{v} = schm0(:,:,s);
        end
        mdl{D} = cat(1,mdl{D},...
            script_inferPH(dt,states,expT,dt_bin,T,repmat(D,1,V),schm,''));
        
        LogL_v = [];
        for v = 1:V
            LogL_v = cat(1,LogL_v,mdl{D}.logL(v));
        end
        logL{D}(:,s) = LogL_v;
    end
end

% calculate BIC for each DPH
Nv = mdl{1}(1).N;
BIC = cell(1,Dmax);
BIC{D} = Inf(V,S);
for D = 1:Dmax
    S = numel(mdl{D});
    for s = 1:S
        df = (D-1)+sum(sum(mdl{D}(s).schm{1})); % (D-1) initial prob, trans prob
        for v = 1:V
            BIC{D}(v,s) = df*log(Nv(v))-2*logL{D}(v,s);
        end
    end
end

% select minimum BIC
degen = zeros(1,V);
minBIC = repmat([Inf,0,0],V,1);
for D = 1:Dmax
    for v = 1:V
        [val,s] = min(BIC{D}(v,:));
        if val<minBIC(v,1)
            minBIC(v,:) = [val,D,s];
        end
    end
end

% show most sufficient state configuration
id = [];
for v = 1:V
    degen(v) = minBIC(v,2);
    id = cat(2,id,repmat(v,[1,degen(v)]));
end
fprintf(['Most sufficient state configuration:\n[%0.2f',...
    repmat(',%0.2f',[1,numel(states(id))-1]),']\n'],states(id));

% infer "true" DPH parameters
disp('Optimize DPH distributions on authentic dwell time histograms...')
schm = cell(1,V);
for v = 1:V
    schm{v} = mdl{minBIC(v,2)}(minBIC(v,3)).schm{v};
end
mdl = script_inferPH(dt,states,expT,1,T,degen,schm,savecurve);

% save computation time
mdl.t_dphtest = toc(t_comp);

fprintf('Most sufficient model complexity found in %0.0f seconds\n',...
    toc(t_comp));


function schm = addExitProb2TransSchemes(schm)

J = size(schm,1);
for nb1 = 1:J
    for id1 = 1:nb1
        for col = 1:J
        end
    end
end


