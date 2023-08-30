function DPHtest_analysisRoutine(rootdir,destdir,varargin)
% DPHtest_analysisRoutine(rootdir,destdir)
% DPHtest_analysisRoutine(rootdir,destdir,subfolders)
%
% Perform simulation and analysis routine used in the main article.
%
% rootdir: source directory where data folder are present
% destdir: destination directory where the results summary is exported
% subfolders: data sub-folders to be analyzed

% defaults
R = 3; % number of simulation replicates
rate = 10; % simulated frame rate
dphprm.excl = false;
dphprm.T = 5;
T_bw = 5;

% set MATLAB search path
codePath = fileparts(mfilename('fullpath'));
addpath(genpath(codePath));

% collect subfolders
if ~isempty(varargin)
    subfolders = varargin{1};
    issubdir = true;
else
    subfolders = {'dataset1','dataset2','dataset3','dataset4','dataset5',...
        'dataset6','dataset7','dataset8','dataset9','dataset10',...
        'dataset11','dataset12','EBS-IBS'};
    issubdir = false;
end

% compile mex files
checkMASHmex();

% create simulation presets files
if issubdir
    DPHtest_createSimPrm(rootdir,subfolders);
else
    DPHtest_createSimPrm(rootdir);
end

% loop analysis through root folder's sub-folders
dircnt = dir(rootdir);
for r = 1:R
    for d = 1:size(dircnt,1)
        if any(contains({'.','..'},dircnt(d,1).name)) || ...
                (issubdir && ~any(cellstrcmp(subfolders,dircnt(d,1).name)))
            continue
        end
        analyzeDirContent(dircnt(d,1),rate,r,R,dphprm,T_bw);
    end
end

% collect results
DPHtest_collectResults(rootdir,destdir);


function analyzeDirContent(cnt,rate,r,R,dphprm,T)

% default
direbsibs = 'trajectories';

% filter out non-data folders
if cnt.isdir && ~strcmp(cnt.name,direbsibs)
    cnt = dir([cnt.folder,filesep,cnt.name]);
    for d = 1:size(cnt,1)
        if any(contains({'.','..'},cnt(d,1).name))
            continue
        end
        analyzeDirContent(cnt(d,1),rate,r,R,dphprm,T);
    end
    
else
    % select _simprm.mat preset files
    issim = true;
    if cnt.isdir && strcmp(cnt.name,direbsibs) && r==1
        setname = 'EBS-IBS';
        subdir = cnt.folder;
        issim = false;
    elseif ~(endsWith(cnt.name,'_simprm.mat'))
        return
    else
        [~,fname,~] = fileparts(cnt.name);
        setname = fname(1:end-length('_simprm'));
        subdir = [cnt.folder,filesep,setname];
    end
    
    % determine whether BW analysis must be performed
    bwana = any(cellstrcmp(setname,{'presets_quench','presets_dyn',...
        'presets_connex','EBS-IBS'}));
    
    % data simulation/import and export
    if issim
        simresfle = [subdir,filesep,setname,sprintf('_%i_simres.mat',r)];
        if any(cellstrcmp(setname,...
                {'presets_quench','presets_dyn','presets_connex'}))
            simtraj = true;
        else
            simtraj = false;
        end
        if exist(simresfle,'file')
            disp(['read simulation parameters from file: ',setname,...
                '_simprm.mat ...']);
            simprm = load([cnt.folder,filesep,setname,'_simprm.mat']);

            disp(['read simulated data from files:',setname,...
                sprintf('_%i_simres.mat',r),' ...']);
            res = load(simresfle);
            res = res.res;
            if simtraj && ~isfield(res,'seq')
                N = numel(res.dt_obs);
                res.seq = cell(1,N);
                for n = 1:N
                    res.seq{n} = ...
                        getDiscrFromDt(res.dt_obs{n},1/simprm.rate);
                end
                save(simresfle,'res','-mat');
            end
            simres = res;
        else
            disp(['simulate data for dataset: ',setname,' ...']);
            [simprm,simres] = DPHtest_simdata(...
                [cnt.folder,filesep,cnt.name],simresfle,simtraj);
        end
    else
        expdatfle = [cnt.folder,filesep,setname,'_data.mat'];
        [simprm,simres] = DPHtest_importexpdata(expdatfle,...
            [cnt.folder,filesep,cnt.name]);
    end
    if isempty(simres)
        return
    end
    
    % data analysis/import and export
    [~,datadir,~] = fileparts(cnt.folder);
    if any(cellstrcmp({'dataset1','dataset2','dataset3','dataset10',...
            'dataset11'},datadir))
        dphprm.Dmax = 4;
        dphprm.bin = 1;
    elseif strcmp('dataset12',datadir)
        dphprm.Dmax = 5;
        dphprm.bin = 1;
    elseif strcmp(datadir,'EBS-IBS')
        dphprm.Dmax = 3;
        dphprm.bin = 10;
    else
        dphprm.Dmax = 3;
        dphprm.bin = 1;
    end

    % import exisiting results
    mldphresfle = [subdir,filesep,setname,'_',num2str(r),'_mldphres.mat'];
    if exist(mldphresfle,'file') && bwana
        mldphres = load(mldphresfle);
        if isfield(mldphres,'fulldphres')
            fulldphres = mldphres.fulldphres;
        end
    end

    % ML-DPH analysis
    if ~exist(mldphresfle,'file') || ...
            (bwana && ~exist('fulldphres','var'))
        [dphres,expres,fulldphres] = MLDPH_analysis(...
            [subdir,filesep,setname,'_',num2str(r)],simres.dt_obs,...
            simprm,dphprm,false);
        if isempty(dphres) || isempty(expres) || isempty(fulldphres)
            return
        end

        % export DPH fit parameters and computation time
        save(mldphresfle,'dphres','expres','fulldphres','-mat');
    end

    % BW inference on simulated state sequences
    if bwana
        if issim && ~iscorrectdegeneracy(simprm,fulldphres)
            return
        end
        
        bwres_fle = ...
            [subdir,filesep,setname,'_',num2str(r),'_bwres.mat'];
        if exist(bwres_fle,'file')
            return
        end

        % using ML-DPH outcome
        if any(cellstrcmp(setname,{'presets_quench','presets_dyn'}))
            bwres_w = BW_analysis(simprm,dphprm,simres.dt_obs,simres.seq,...
                fulldphres,T,true);
            if isempty(bwres_w)
                return
            end
        else
            bwres_w = [];
        end

        % not using ML-DPH outcome
        bwres_wo = BW_analysis(simprm,dphprm,simres.dt_obs,simres.seq,...
            fulldphres,T,false);
        if isempty(bwres_wo)
            return
        end

        % save results to _bwres.mat file
        save(bwres_fle,'bwres_w','bwres_wo','-mat');
    end
end


function ok = iscorrectdegeneracy(simprm,dphres)

val0 = unique(sort(simprm.val));
V0 = numel(val0);
D0 = zeros(1,V0);
for v = 1:V0
    D0(v) = sum(simprm.val==val0(v));
end

val = unique(sort(dphres{4}));
V = numel(val);
D = zeros(1,V);
for v = 1:V
    D(v) = sum(dphres{4}==val(v));
end

ok = all(D0==D);
