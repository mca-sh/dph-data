function [prm,res] = PHtest_simdata(flepresets,fleout,simtraj,PHtype,isiamm)
% [prm,res] = PHtest_simdata(filepresets,simtraj,PHtype)
%
% Simulate replicate data sets from input presets file.
% 
% filepresets: path to simulation parameter file (*_simprm.mat)
% fleout: destination file (*_simres.mat)
% simtraj: (1) simulate state sequences, (0) simulate only dwell times
% PHtype: distribution type if simulated dwell times are drawn directly 
%         from distribution 1 for discrete PH, 2 for continuous PH
% prm: simulation parameters
% res: simulated dwell times and trajectories
%
% example: 
% >> fleprm = 'dataset2/presets_II_13_simprm.mat';
% >> fleres = 'dataset2/presets_II_13/presets_II_13_1_simres.mat';
% >> [prm,res] = PHtest_simdata(fleprm,fleres,false,2)

% defaults (from iAMM.m Hines 2015)
E_mu = [0,0.7];      % state means 
E_sig2 = [0.01; 0.01]; % state variances 

% collect simulation parameters
prm = load(flepresets);

% create destination folder if none
[dest,~,~] = fileparts(fleout);
if ~exist(dest,'dir')
    mkdir(dest);
end

% generate data
if isiamm
    stateval = unique(prm.val);
    nAgg = numel(stateval);
    J = numel(prm.val);
    E.mu = zeros(1,J);
    E.sigma2 = zeros(1,J);
    obs_id = zeros(1,J);
    for agg = 1:nAgg
        E.mu(prm.val==stateval(agg)) = E_mu(agg);
        E.sigma2(prm.val==stateval(agg)) = E_sig2(agg);
        obs_id(prm.val==stateval(agg)) = agg;
    end
    A = prm.tp;
    A(~~eye(J)) = 1-sum(prm.tp,2);
    res.traj = cell(1,prm.N);
    res.seq = cell(1,prm.N);
    res.seq_gt = cell(1,prm.N);
    res.dt_obs = cell(1,prm.N);
    res.dt_gt = cell(1,prm.N);
    for n = 1:prm.N
        [Y_n,S_GT_n] = HmmGenerateData_MCASH(prm.ndt, prm.ip', A, E, 'normal');
        res.traj{n} = Y_n;
        res.seq{n} = obs_id(S_GT_n);
        res.seq_gt{n} = S_GT_n;

        res.dt_obs{n} = getDtFromDiscr(res.seq{n}',1);
        res.dt_obs{n}(end,3) = res.dt_obs{n}(end,2); % remove NaN
        res.dt_obs{n} = [res.dt_obs{n},stateval(res.dt_obs{n}(:,[2,3]))];
        res.dt_obs{n}(:,1) = res.dt_obs{n}(:,1)/prm.rate;

        res.dt_gt{n} = getDtFromDiscr(res.seq_gt{n}',1);
        res.dt_gt{n}(end,3) = res.dt_gt{n}(end,2); % remove NaN
        res.dt_gt{n} = [res.dt_gt{n},E.mu(res.dt_gt{n}(:,[2,3]))];
        res.dt_gt{n}(:,1) = res.dt_gt{n}(:,1)/prm.rate;
    end
    res.nAgg = nAgg;
else
    res = PHtest_buildModel(prm,simtraj,PHtype);
end

% export simulated dwell times to mat file
save(fleout,'res','-mat');

