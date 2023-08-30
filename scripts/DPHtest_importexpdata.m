function [prm,dat] = DPHtest_importexpdata(destfle,src)
% [prm,dat] = DPHtest_importData(rootdir)
%
% Import state trajectories and dwell times from .traces files.
%
% src: source directory
% destfle: path to destination data file

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
N = size(flist,1);

% read trajectories
dat.dt_obs = cell(1,N);
dat.seq = cell(1,N);
states0 = [];
for n = 1:N
    % import time and FRET state data
    traj = importdata([src,flist(n,1).name],'\t',3);
    dat.seq{n} = traj.data(:,end);
    
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
V = numel(prm.val);

% add state indexes to dwell time table
for n = 1:N
    stateid = zeros(size(dat.dt_obs{n},1),2);
    for v = 1:V
        stateid(dat.dt_obs{n}(:,[2,3])==prm.val(v)) = v;
    end
    stateid(stateid==0) = NaN;
    dat.dt_obs{n} = ...
        cat(2,dat.dt_obs{n}(:,1),stateid,dat.dt_obs{n}(:,2:end));
    for v = 1:V
        dat.seq{n}(dat.seq{n}==prm.val(v)) = v;
    end
end

% save to file
save(destfle,'prm','dat','-mat');
