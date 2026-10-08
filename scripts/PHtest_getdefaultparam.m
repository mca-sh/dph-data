function prm = PHtest_getdefaultparam()
% prm = PHtest_getdefaultparam()
%
% Defines default parameters specific to ML-DPH analysis and return them in
% a structure.
%
% prm: structure containing analysis parameters with fields:
%   prm.T_BW: nb. of ML-BW restarts
%   prm.excl: 1 to exclude first & last dwell times of each trajectory, 0
%             otherwise.
%   prm.bin: bin size
%   prm.Dmin: min. number of degenerate states to fit
%   prm.Dmax: max. number of degenerate states to fit
%   prm.applyrules: apply additional selection rules on inferred model 
%                   (distribution validity and state doublons)

prm.T_BW = 5;           % number of ML-BW restarts
prm.bin = 1;            % dwell-time binning (in # of time steps)
prm.excl = false;       % exclude first and last dwell times from histogram
prm.applyrules = false; % apply additional selection rules on inferred model
prm.Dmin = [1,1];       % minimum nb of states tested
prm.Dmax = [4,4];       % maximum nb of states tested