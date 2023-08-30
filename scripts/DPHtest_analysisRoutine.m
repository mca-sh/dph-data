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
dphprm.bin = 1;
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
        'dataset6','dataset7'};
    issubdir = false;
end

% compile mex files
checkMASHmex();

% create simulation presets files
disp('Generate simulation presets files...');
if issubdir
    DPHtest_createSimPrm(rootdir,subfolders);
else
    DPHtest_createSimPrm(rootdir);
end

% loop analysis through root folder's sub-folders
dircnt = dir(rootdir);
for d = 1:size(dircnt,1)
    if any(contains({'.','..'},dircnt(d,1).name)) || ...
            (issubdir && ~any(contains(subfolders,dircnt(d,1).name)))
        continue
    end
    analyzeDirContent(dircnt(d,1),rate,R,dphprm,T_bw);
end

% collect results
DPHtest_collectResults(rootdir,destdir);


function analyzeDirContent(cnt,rate,R,dphprm,T)

% filter out non-data folders
if cnt.isdir
    cnt = dir([cnt.folder,filesep,cnt.name]);
    for d = 1:size(cnt,1)
        if any(contains({'.','..'},cnt(d,1).name))
            continue
        end
        analyzeDirContent(cnt(d,1),rate,R,dphprm,T);
    end
    
else
    % select _simprm.mat preset files
    if ~(endsWith(cnt.name,'_simprm.mat'))
        return
    end
    [~,fname,~] = fileparts(cnt.name);
    setname = fname(1:end-length('_simprm'));
    subdir = [cnt.folder,filesep,setname];
    
    % determine whether BW analysis must be performed
    bwana = any(contains({'presets_quench','presets_dyn'},setname));
    
    % data simulation/import and export
    simdatfles = cell(1,R);
    str_fles = '';
    for r = 1:R
        simdatfles{r} = ...
            [subdir,filesep,setname,sprintf('_%i_simres.mat',r)];
        str_fles = cat(2,str_fles,' ',setname,sprintf('_%i_simres.mat',r));
    end
    if  ~any(cellfun(@isempty,simdatfles)) && ...
            all(cellfun(@exist,simdatfles))
        disp(['read simulation parameters from file: ',setname,...
            '_simprm.mat ...']);
        simprm = load([cnt.folder,filesep,setname,'_simprm.mat']);
        
        disp(['read simulated data from files:',str_fles,' ...']);
        simres = cell(1,R);
        for r = 1:R
            res = load(simdatfles{r});
            simres{r} = res.res;
        end
    else
        disp(['simulate data for dataset: ',setname,' ...']);
        if any(contains({'presets_quench','presets_dyn'},setname))
            simtraj = true;
        else
            simtraj = false;
        end
        [simprm,simres] = ...
            DPHtest_simdata([cnt.folder,filesep,cnt.name],R,simtraj);
    end
    
    % data analysis/import and export
    [~,datadir,~] = fileparts(cnt.folder);
    if any(contains({'dataset1','dataset2'},datadir))
        dphprm.Dmax = 5;
    else
        dphprm.Dmax = 4;
    end
    for r = 1:R
        if isempty(simres{r})
            continue
        end

        mldphres_fle = ...
            [subdir,filesep,setname,'_',num2str(r),'_mldphres.mat'];
        
        % import exisiting results
        if exist(mldphres_fle,'file') && bwana
            mldphres = load(mldphres_fle);
            dphres = mldphres.dphres;
            
        % ML-DPH analysis
        elseif ~exist(mldphres_fle,'file')
            dphres = MLDPH_analysis(...
                [subdir,filesep,setname,'_',num2str(r),'_dph'],...
                simres{r}.dt_obs,simprm,dphprm,false);
            if isempty(dphres)
                continue
            end

            expres = MLDPH_analysis(...
                [subdir,filesep,setname,'_',num2str(r),'_exp'],...
                simres{r}.dt_obs,simprm,dphprm,true);
            if isempty(expres)
                continue
            end

            % export DPH fit parameters and computation time
            save(mldphres_fle,'dphres','expres','-mat');
        end
        
        % BW inference on simulated state sequences
        if bwana
            bwres_fle = ...
                [subdir,filesep,setname,'_',num2str(r),'_bwres.mat'];
            if exist(bwres_fle,'file')
                continue
            end
            
            % using ML-DPH outcome
            bwres_w = BW_analysis(simprm,dphprm,simres{r}.dt_obs,...
                simres{r}.seq,dphres,T,true);
            if isempty(bwres_w)
                continue
            end

            % not using ML-DPH outcome
            bwres_wo = BW_analysis(simprm,dphprm,simres{r}.dt_obs,...
                simres{r}.seq,dphres,T,false);
            if isempty(bwres_wo)
                continue
            end

            % save results to _bwres.mat file
            save(bwres_fle,'bwres_w','bwres_wo','-mat');
        end
    end
end
