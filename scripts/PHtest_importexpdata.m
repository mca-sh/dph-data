function [prm,dat] = PHtest_importexpdata(destfle,src)
% [prm,dat] = PHtest_importData(rootdir)
%
% Import state trajectories and dwell times from .traces files.
%
% src: source directory
% destfle: path to destination data 

% CONSTANTS
THRESH = 0.5; % threshold between real states
STATEVAL = [0,0.7]; % state values to associate for .txt traj files (newly discetized with vbFRET)
BS_NDTMIN = 0; % min. nb of dwell times after bootstrapping (0 = no BS)

% check if data were already imported
if exist(destfle,'file')
    flecnt = load(destfle);
    prm = flecnt.prm;
    dat = flecnt.dat;
    return
end

% prop'ups path
if src(end)~=filesep
    src = [src,filesep];
end

% list traces files
flist = dir([src,'*.traces']);
if isempty(flist)
    flist = dir([src,'*.txt']);
end
N = size(flist,1);

% read trajectories
dat.traj = cell(1,N);
dat.dt_obs = cell(1,N);
dat.seq = cell(1,N);
states0 = [];
for n = 1:N
    % import FRET and state data
    traj = importdata([src,flist(n,1).name],'\t',3);
    dat.traj{n} = traj.data(:,end-1);
    dat.seq{n} = traj.data(:,end);
    
    % associate states to reference values for .txt files
    [~,~,fext] = fileparts(flist(n,1).name);
    if strcmp(fext,'.txt')
        dat.seq{n}(dat.seq{n}<THRESH) = STATEVAL(1);
        dat.seq{n}(dat.seq{n}>=THRESH) = STATEVAL(2);
    end
    
    % calculates sampling frame rate
    if n==1
        prm.rate = 1/(traj.data(2,1)-traj.data(1,1));
    end
    
    % collect state values
    states0 = cat(2,states0,unique(dat.seq{n}'));
    
    % calculates dwell times
    dat.dt_obs{n} = getDtFromDiscr(dat.seq{n},1/prm.rate);
end
prm.val = sort(unique(states0));

% add state indexes to dwell time table and count dwell times for each 
% aggregate
V = numel(prm.val);
ndt = zeros(N,V);
for n = 1:N
    stateid = zeros(size(dat.dt_obs{n},1),2);
    for v = 1:V
        stateid(dat.dt_obs{n}(:,[2,3])==prm.val(v)) = v;
        ndt(n,v) = nnz(dat.dt_obs{n}(:,2)==prm.val(v));
    end
    stateid(stateid==0) = NaN;
    dat.dt_obs{n} = ...
        cat(2,dat.dt_obs{n}(:,1),stateid,dat.dt_obs{n}(:,2:end));
    for v = 1:V
        dat.seq{n}(dat.seq{n}==prm.val(v)) = v;
    end
end

% bootstrap trajectrories to artificially construct a sample with suffiient 
% stats
while any(sum(ndt,1)<BS_NDTMIN)
    % draw randomly a molecule from the dataset
    n = randi(N,1);
    
    % append the original sample with molecule's data
    dat.traj = cat(2,dat.traj,dat.traj{n});
    dat.dt_obs = cat(2,dat.dt_obs,dat.dt_obs{n});
    dat.seq = cat(2,dat.seq,dat.seq{n});
    ndt = cat(1,ndt,ndt(n,:));
end

dat.nAgg = V;

% save to file
save(destfle,'prm','dat','-mat');
