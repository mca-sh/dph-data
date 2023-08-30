function [prm,res] = DPHtest_simdata(flepresets,fleout,simtraj)
% [prm,res] = DPHtest_simdata(filepresets,simtraj)
%
% Simulate replicate data sets from input presets file.
% 
% filepresets: path to presets file (*.mat)
% fleout: destination file
% simtraj: (1) simulate state sequences, (0) simulate only dwell times
% prm: simulation parameters
% res: simulated dwell times and trajectories
%
% example: 
% >> [prm,res] = DPHtest_simdata(dataset2/presets_II_13.mat,10,3)

% collect simulation parameters
prm = load(flepresets);

% create destination folder if none
[dest,~,~] = fileparts(fleout);
if ~exist(dest,'dir')
    mkdir(dest);
end

% generate data
res = DPHtest_buildModel(prm,simtraj);

% export simulated dwell times to mat file
save(fleout,'res','-mat');

