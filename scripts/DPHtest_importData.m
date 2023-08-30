function [prm,dat] = DPHtest_importData(rootdir)
% [prm,dat] = DPHtest_importData(rootdir)
%
% Import state trajectories and dwell times from .traces files.
%
% rootdir: source directory

prm = [];
dat = [];

if ~exist(rootdir,'dir')
    disp('DPHtest_importData: source directory not found.');
end
if rootdir(end)~=filesep
    rootdir = [rootdir,filesep];
end
flst = dir([rootdir,'*.traces']);
for f = 1:size(flst,1)
    flst(f,1).name
end