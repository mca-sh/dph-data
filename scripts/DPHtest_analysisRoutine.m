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
    
    % select .mat preset files
    if ~(endsWith(cnt.name,'_simprm.mat'))
        return
    end
    [~,fname,~] = fileparts(cnt.name);
    setname = fname(1:end-length('_simprm'));
    subdir = [cnt.folder,filesep,setname];
    
    % show action
    disp(['simulate and analyze dataset: ',setname]);
    
    % data simulation and export
    if any(contains({'presets_quench','presets_dyn'},setname))
        simtraj = true;
    else
        simtraj = false;
    end
    [simprm,simres] = ...
        DPHtest_simdata([cnt.folder,filesep,cnt.name],R,simtraj);
    
    % data analysis and export
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

        % ML-DPH analysis of simulated dwell times
        dphres = ...
            MLDPH_analysis([subdir,filesep,setname,'_',num2str(r),'_dph'],...
            simres{r}.dt_obs,simprm,dphprm,false);
        if isempty(dphres)
            continue
        end
        
        expres = ...
            MLDPH_analysis([subdir,filesep,setname,'_',num2str(r),'_exp'],...
            simres{r}.dt_obs,simprm,dphprm,true);
        if isempty(expres)
            continue
        end

        % export DPH fit parameters and computation time
        fname = [setname,'_',num2str(r),'_mldphres.mat'];
        save([subdir,filesep,fname],'dphres','expres','-mat');
        
        % BW inference on simulated state sequences
        if any(contains({'presets_quench','presets_dyn'},setname))
            
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
            fname = [setname,'_',num2str(r),'_bwres.mat'];
            save([subdir,filesep,fname],'bwres_w','bwres_wo','-mat');
        end
    end
end
