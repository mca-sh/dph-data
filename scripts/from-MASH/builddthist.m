function [dt,prob,edg] = builddthist(dttbl,val,rate,bin,excl)
% [dt,prob,edg] = builddthist(dttbl,val,rate,bin,excl)
%
% Build dwell time histogram from single molecule dwell time tables as in 
% MLPH_analysis and script_findBestModel.
%
% dttbl: {1-by-N}[ndt_n-by-5] dwell time, states, state values
% val: [1-by-J] FRET values
% rate: trajectory sampling rate
% bin: dwell time binning
% excl: (1) exclude first and last dwell times of state sequences, (0)
%       otherwise
% dt: [1-by-nbin] histogram bins
% prob: [1-by-nbin]: histogram counts normalized by the sum
% edg: [1-by-(nbin+1)] bin edges

% get state sequences
N = numel(dttbl);
dt_new = [];
for n = 1:N
    dt_m = dttbl{n};
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
dttbl = dt_new;
if isempty(dttbl)
    disp(['ML-PH can not proceed: no dwell times are left after ',...
        'exclusion of static trajectories.']);
    return
end

valv = unique(val);
V = numel(valv);
dt = cell(1,V);
prob = cell(1,V);
edg = cell(1,V);
for v = 1:V
    dt_v = dttbl(dttbl(:,2)==v,:);
    if isempty(dt_v)
        continue
    end
    dt_v(:,1) = dt_v(:,1)*rate;
    
    % check if dt table contains original PDF
    ispdf = size(dttbl,2)==4; % time, state1, 0, prob
    
    % binned dwell time histogram
    if bin~=1
        edg{v} = 1:bin:(max(dt_v(:,1))+bin);
        x = floor(mean([edg{v}(2:end)-1;edg{v}(1:end-1)],1));
    else
        edg{v} = 1:(max(dt_v(:,1))+1); % [edg1,edg2[
        x = edg{v}(1:end-1);
    end
    [P,~,id] = histcounts(dt_v(:,1),edg{v});
    if ispdf 
        P = zeros(size(P));
        for i = 1:numel(P)
            P(i) = sum(dt_v(id==i,4));
        end
    end
    P = P/sum(P);
%     validprob = P>0;
%     dt{v} = x(validprob);
%     prob{v} = P(validprob);
%     edg{v}(find(~validprob)+1) = [];
    dt{v} = x;
    prob{v} = P;
end

