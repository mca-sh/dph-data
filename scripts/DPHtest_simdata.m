function [prm,allres] = DPHtest_simdata(filepresets,R,simtraj)
% [prm,res] = DPHtest_simdata(filepresets,R,simtraj)
%
% Simulate replicate data sets from input presets file.
% 
% filepresets: path to presets file (*.mat)
% R: number of replicates
% simtraj: (1) simulate state sequences, (0) simulate only dwell times
% prm: simulation parameters
% res: simulated dwell times and trajectories
%
% example: 
% >> [prm,res] = DPHtest_simdata(dataset2/presets_II_13.mat,10,3)

% defaults
allres = cell(1,R);

% collect simulation parameters
prm = load(filepresets);

% get data set name
[pname,fname,~] = fileparts(filepresets);
fname = fname(1:end-length('_simprm'));

% generate and export data nRep times
if pname(end)~=filesep
    pname = [pname,filesep];
end
pname_out = [pname,fname];
if ~exist(pname_out,'dir')
    mkdir(pname_out);
end
if pname_out(end)~=filesep
    pname_out = [pname_out,filesep];
end
for r = 1:R
    if R>1
        fname_out = [fname,'_',num2str(r),'_simres'];
    else
        fname_out = [fname,'_simres'];
    end

    % generate data
    res = DPHtest_buildModel(prm,simtraj);
    allres{r} = res;

    % export simulated dwell times to mat file
    save([pname_out,fname_out,'.mat'],'res','-mat');
end


