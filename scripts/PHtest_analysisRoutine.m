function PHtest_analysisRoutine(src0,figdir,meth,varargin)
% PHtest_analysisRoutine(src, dest, method_code)
% PHtest_analysisRoutine(_, specific_data_folders)
% PHtest_analysisRoutine(_, dump_dir_name)
%
% Perform simulation and ML-DPH, iEMM and ML-ENMM analysis routine used in 
% article "Using Phase-type Distributions to Model Kinetic Heterogeneity in 
% Ensemble Dwell time histograms", Hadzic et. el. 2026
%
% src         : source directory where data folder are present
% dest        : figure directory where the result summaries are exported
% method_code :-2: EM-EMM (with positive coefficients)
%              -1: EM-ENMM (with negative coefficients)w
%               0: test canonical connectivities (+10 for searching for 
%                  min. complexity)
%               1: test all connectivities (+10 for searching for min. 
%                  complexity)
%               2: test "uncoupled" and "irreversible loop" 
%                  connectivities (+10 for searching for min. complexity)
%               3: test "uncoupled" connectivity (+10 for searching for 
%                  min. complexity)
%               4: test "coupled"  (+10 for searching for min. 
%                  complexity)
%               5: test "uncoupled" and "generalized coxian" connectivity 
%                  (+10 for searching for min. complexity)
%               6: test "acyclic" connectivity  (+10 for searching for 
%                  min. complexity).
%               7: test "acyclic" connectivity with different state 
%                  initiations.
%              10: iEMM adapted from Hines et. al. 2015
% specific_data_folders: {1-by-nFolders} data sub-folders to be analyzed
% dump_dir_name: char dump directury name

    % CONSTANTS
    R = 10; % number of simulation replicates
    
    % set MATLAB search path
    codePath = fileparts(mfilename('fullpath'));
    addpath(genpath(codePath));
    
    % propup main paths
    if src0(end)~=filesep
        src0 = [src0,filesep];
    end
    if figdir(end)~=filesep
        figdir = [figdir,filesep];
    end
    
    % compile mex files
    checkMASHmex();
    
    % list main "dataset" folders located in source directory
    [d_0, specfldr, isspec, dumpname] = PHtest_list_data_folders(src0, ...
        varargin);
    
    diary off % stop logging to not overload .log file
    
    % simulate or import data and list all input/output files
    [f_prm, datnm, f_dat_r, f_mlph, f_bw, dat_id, phprm0] = ...
        PHtest_initialization(src0, meth, R, d_0, specfldr, isspec, ...
        dumpname);
    
    % Initiate printing (necessary for a large list of files)
    disp('Analyze datasets ...');
    [n_exist, n_all, dat_name, max_l, ndigit] = init_print_progress('', ...
        datnm, f_mlph, dat_id);
    
    % Print the initial list nb of completed analysis
    nchar = print_progress_parallel(dat_name, n_exist, n_all, max_l, ...
        ndigit, 0);
    
    % Define data queue function for progress printing during parallelized
    % process
    q = parallel.pool.DataQueue();
    afterEach(q, @(id) parfor_loop_callback(id));
        function parfor_loop_callback(id)
            % Increment on the main thread the nb of completed analysis for 
            % the input dataset
            n_exist(id) = n_exist(id) + 1;

            % Print current state
            print_progress_parallel(dat_name, n_exist, n_all, max_l, ...
                ndigit, nchar);
        end
    
    % Main analysis parallelized loop
    N = size(f_dat_r,1);
    parfor n = 1:N 
    
        % Reset print refresh
        refresh_print = false;
    
        % Get dataset name
        f_dat_n = f_dat_r{n};
        [src_n,~] = fileparts(f_dat_n);
        folder_names = split(src_n((length(src0)+1):end),filesep);
        datdir_n = folder_names{2};
    
        % adjust analysis parameters for current dataset
        phprm = PHtest_adjustparam(phprm0,datdir_n);

        if endsWith(f_dat_n,'_simres.mat')
            % Imports simulated data form file
            sim = load(f_dat_n,'res');
            if ~isfield(sim,'res')
                continue
            end
            dat = sim.res;
            simprm = load(f_prm{n});
        else
            % Imports experimental data form file
            [simprm,dat] = PHtest_importexpdata(f_dat_n,'');
        end
        
        if ~exist(f_mlph{n},'file')
            % Enable degeneracy analysis (ML-DPH, ML-ENNM or iEMM) if
            % results file does not exist
            runana = true;
        elseif contains('phprm', who('-file',f_mlph{n}))
            res = load(f_mlph{n},'dphres','phprm'); % load analysis options "phprm"
            if ~isequal(res.phprm, phprm)
                % Enable degeneracy analysis if analysis options have 
                % changed since last saving
                runana = true;
            else
                % Skip degeneracy analysis if already performed and saved
                dphres = res.dphres;
                runana = false;
            end
        else
            runana = true;
        end

        if runana
            % Permform ML-DPH, ML-ENMM or iEMM analysis on dwell times
            dphres = PHtest_MLPHanalysis('', dat.dt_obs, simprm, phprm, ...
                false);
        
            % Save dwell time analysis results to MATLAB file
            dat2save = struct('phprm',phprm);
            dat2save.dphres = dphres;
            save(f_mlph{n},'-mat','-fromstruct',dat2save);
    
            refresh_print = true;
        end
    
        if ~isempty(f_bw{n})
            if exist(f_bw{n},'file')
                % Skip BW analysis if already performed and saved to file
                bwana = false;
            elseif isempty(f_prm{n})
                % Enable BW analysis for experimental datasets
                bwana = true;
            else
                % Enable BW analysis for simulated datasets whose ground 
                % truth model complexity was correctly retrieved
                bwana = iscorrectdegeneracy(simprm,dphres);
            end

            if bwana
                % Run BW on trajectories using optimum model complexity
                bwres_w = PHtest_BWanalysis(simprm,phprm,dat.dt_obs,dat.seq,...
                    dphres,phprm.T_BW,false,false);
        
                % Save BW results to MATLAB file
                dat2save = struct();
                dat2save.bwres_w = bwres_w;
                save(f_bw{n},'-mat','-fromstruct',dat2save);
            end
        end
    
        if refresh_print
            % Show progress using data queue function
            send(q, dat_id(n)); % Send ID to main thread
        end
    end

    % Resume logging
    diary on 

    % Identify which figure data mist be updated
    fig_labels = fig_labels_from_dataset_names(specfldr);
    
    % Reshape results in a format readible by figure-building scripts
    % (e. g. by DPHartcl_figure3, DPHartcl_figure4_v2, etc.)
    fprintf('Collecting results...\n') 
    PHtest_collectResults(src0, figdir, dumpname, fig_labels, false);
    fprintf('Process completed!\n');
end


function ok = iscorrectdegeneracy(simprm,phres)
% ok = iscorrectdegeneracy(simprm,phres)
%
% Determine whether the GT state degeneracy was recovered by ML-PH.
%
% simprm: structure containing simulation parameters with field:
%   simprm.val: [1-by-J] FRET state values in GT configuration
% phres: {1-by-4} ML-PH results with:
%   phres{4}: [1-by-J'] FRET state values in inferred configuration
% ok: 1 if GT degeneracy was recovered, 0 otherwise
    
    val0 = unique(sort(simprm.val));
    V0 = numel(val0);
    D0 = zeros(1,V0);
    for v = 1:V0
        D0(v) = sum(simprm.val==val0(v));
    end
    
    D = determinestatedegeneracy(phres);
    
    ok = all(D0==D);
end


function D = determinestatedegeneracy(phres)
    val = unique(sort(phres{4}));
    V = numel(val);
    D = zeros(1,V);
    for v = 1:V
        D(v) = sum(phres{4}==val(v));
    end
end


function fig_labels = fig_labels_from_dataset_names(specfldr)
    FIG_LABEL_ALL = {'PERFONSIM1','PERFONSIM2','EBSIBS','SHAPESUCCESS'};
    SET_NAMES_ALL = {...
    {'dataset200', 'dataset201', 'dataset202', 'dataset210', 'dataset211', ...
     'dataset212'} ...
    {'dataset100', 'dataset101', 'dataset102', 'dataset200', 'dataset201', ...
     'dataset202', 'dataset300', 'dataset301', 'dataset302', 'dataset400', ...
     'dataset401', 'dataset402'} ...
    {'EBS-IBS'}, ... 
    {'dataset200', 'dataset201', 'dataset202', 'dataset300', 'dataset301', ...
    'dataset302', 'dataset400', 'dataset401', 'dataset402'}...
    };

    num_labels = length(FIG_LABEL_ALL);
    incl = false(1,num_labels);
    for fig = 1:num_labels
        if any(contains(SET_NAMES_ALL{fig},specfldr))
            incl(fig) = true;
        end
    end
    fig_labels = FIG_LABEL_ALL(incl);
end
