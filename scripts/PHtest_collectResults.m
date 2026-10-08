function PHtest_collectResults(rootdir,destdir,varargin)
% PHtest_collectResults(rootdir,destdir)
% PHtest_collectResults(...,subdir)
% PHtest_collectResults(...,dumpdir)
% PHtest_collectResults(...,verbose)
%
% Collect results used in article's figures and export them to .txt files.
%
% rootdir: source directory
% destdir: destination directory
% subdir: (opt) {1-by-D} specific subdirectories to export
% dumpdir: specific dump directory
% verbose: true to show logs, false to mute

% CONSTANTS
FIG_LABEL_ALL = {'PERFONSIM1','PERFONSIM2','EBSIBS','SHAPESUCCESS'};
TRAJ_DIR = 'EBS-IBS/trajectories/';
FILE_TRAJ_1 = [TRAJ_DIR, 'IBS-20mM-Mg-3_all185_mol113_post(2).traces'];
FILE_TRAJ_2 = [TRAJ_DIR, 'IBS-20mM-Mg-3_all185_mol113_post.traces'];
SET_NAMES_ALL = {...
    {   {'dataset200', 'dataset201', 'dataset202'}, ... % PERFONSIM1
        {'dataset210', 'dataset211', 'dataset212'} ...
    }, ... 
    {   {'dataset100', 'dataset101', 'dataset102', 'dataset200', ... % PERFONSIM2
         'dataset201', 'dataset202', 'dataset300', 'dataset301', ...
         'dataset302', 'dataset400', 'dataset401', 'dataset402'} ...
    }, ... 
    {   {FILE_TRAJ_1, FILE_TRAJ_2}, ... % EBSIBS
        {'EBS-IBS'} ...
    }, ... 
    {   {'dataset200', 'dataset201', 'dataset202', 'dataset300', ... % SHAPESUCCESS
         'dataset301', 'dataset302', 'dataset400', 'dataset401', ...
         'dataset402'} ...
    }}; 
VAR_NAMES_ALL = { ...
    {   {'schmid', 'schm2', 'hist2', 'PMF2', 'shapeGT2', 'fit2', 'D2', ... % PERFONSIM1
         'tMLPH'}, ...
        {'tau', 'wba', 'PMF2', 'shapeGT2', 'D2', 'spectralGT2'} ...
    }, ... 
    {   {'schmid', 'schm2', 'hist2', 'PMF2', 'fit2', 'shapeGT2', ... % PERFONSIM2
         'shapefit2', 'D2', 'tMLPH'}}, ...
    {   {'traj'}, ... % EBSIBS
        {'hist1', 'fit1', 'BIC1', 'hist2', 'fit2', 'BIC2', 'tMLPH', ...
         'TPMBW', 'states', 'pop', 'tBW'}}, ...
    {   {'schmid', 'schm2', 'PMF2', 'shapeGT2', 'D2', 'spectralGT2'} ... % SHAPESUCCESS
    }}; 

% set MATLAB search path
codePath = fileparts(mfilename('fullpath'));
addpath(genpath(codePath));

% collect input arguments
fig_label = FIG_LABEL_ALL;
dumpname = [];
verbose = true;
for arg = varargin
    if iscell(arg{1})
        fig_label = arg{1};
    elseif ischar(arg{1}) && isempty(dumpname)
        dumpname = arg{1};
    elseif islogical(arg{1})
        verbose = arg{1};
    end
end

% clean input from unknown figure labels
incl = contains(FIG_LABEL_ALL, fig_label);
fig_label = FIG_LABEL_ALL(incl);
set_names_all = SET_NAMES_ALL(incl);
var_names_all = VAR_NAMES_ALL(incl);

% check existence of source directory
if ~exist(rootdir,'dir')
    if verbose
	    disp('Source directory not found.')
    end
	return
end

% propup paths
if rootdir(end)~=filesep
	rootdir = [rootdir,filesep]; % used in eval()
end
if destdir(end)~=filesep
	destdir = [destdir,filesep];
end

% add dump directory
if ~isempty(dumpname)
    destdir = [destdir,'data_',dumpname,filesep];
end

% create destination directory if not existing
if ~exist(destdir,'dir')
    if verbose
        disp('create destination directory...')
    end
    mkdir(destdir);
end

% Collect handles to export functions based on figure labels
num_labels = length(fig_label);
export_functions = cell(1, num_labels);
for sf = 1:num_labels
    export_functions{sf} = str2func(['exportfigdat_', fig_label{sf}]);
end

% show process
if verbose
    strsubfig = '';
    for sd = 1:num_labels
        strsubfig = cat(2,strsubfig,fig_label{sd},' ');
    end
    disp(['files will be exported for: ',strsubfig(1:end-1),'.']);
end

% Gather unique names of variables and datasets to be imported in a 
% preliminary step
set_names = {};
var_names = {};
for sf = 1:numel(fig_label)
    num_lots_sf = length(set_names_all{sf});
    for lot = 1:num_lots_sf
        num_sets_lot = length(set_names_all{sf}{lot});
        var_names_lot = var_names_all{sf}{lot};
        for d = 1:num_sets_lot
            set_name_d = set_names_all{sf}{lot}{d};
            [set_names, var_names] = append_set_and_var_lists(set_names, ...
                var_names, set_name_d, var_names_lot);
        end
    end
end

% Import variables from files for each dataset
num_sets = length(set_names);
dat0 = cell(1,num_sets);
for d = 1:num_sets
    dat0{d} = readData([rootdir, set_names{d}], dumpname, verbose, ...
        var_names{d}{:});
end

% Format and export data in individual figure files
for sf = 1:numel(fig_label)
    if verbose
        disp(['process ',fig_label{sf},'...']);
    end

    num_lots_sf = length(set_names_all{sf});
    dat_sf = cell(1, num_lots_sf);
    for lot = 1:num_lots_sf
        num_sets_lot = length(set_names_all{sf}{lot});
        var_names_lot = var_names_all{sf}{lot};
        dat_sf{lot} = cell(1, num_sets_lot);
        for d = 1:num_sets_lot
            set_name_d = set_names_all{sf}{lot}{d};
            dat_sf{lot}{d} = collect_var_in_list(dat0, set_names, ...
                var_names, set_name_d, var_names_lot);
        end
    end

    export_function = export_functions{sf};
    export_function(dat_sf, destdir, verbose);
end

% Show success
if verbose
    disp('Process completed!');
end


function [set_names, var_names] = append_set_and_var_lists(set_names, ...
    var_names, set_name_d, var_names_lot)
% Append list of datasets and variable names with input if absent from the 
% list.

ds = find(strcmp(set_names, set_name_d), 1);
if isempty(ds)
    set_names = cat(2, set_names, set_name_d);
    var_names = cat(2, var_names, {var_names_lot});
else
    var_id = find(~contains(var_names_lot, var_names{ds}));
    for var = var_id
        var_names{ds} = cat(2, var_names{ds}, var_names_lot{var});
    end
end


function dat = collect_var_in_list(dat0, set_names, var_names, set_name_d, ...
    var_names_d)
% Collect the values of variables with input names for an input dataset
% name.

num_var = length(var_names_d);
dat = cell(1,num_var);

ds = find(strcmp(set_names, set_name_d), 1);
for var = 1:num_var
    dat{var} = dat0{ds}{strcmp(var_names_d{var}, var_names{ds})};
end



function exportfigdat_PERFONSIM1(dat, dest, verbose)
% exportfigdat_PERFONSIM1(dat, dest, verbose)
%
% Restructure and export data for figure PERFONSIM1 of the main article.
%
% dat: {1-by-2} values of variables imported from files
% dest: destination directory
% verbose: true to show logs, false to mute

% CONSTANTS (headers)
H0 = {'distribution_shapes','resolution_power'};
HA1 = {'schm id','schm','hist','PMF_GT','shape_GT','PMF_fit','Db',...
    'computation(s)'};
HA10 = {'dt','normcounts'};
HA11 = {'dt','PMF'};
HA12 = {'dt','fit','cumfit'};
HB1 = {'tau','wba','PMF','shape','Db','spectral_GT'};
HB2 = {'dt','PMF'};

% % show process
% if verbose
%     disp('>> collect data for figure PERFONSIM1...');
% end

% d200 = {get_var_value(DIRA0, set_names, var_names, 'schmid', 'schm2', ...
%     'hist2', 'PMF2', 'shapeGT2', 'fit2', 'D2', 'tMLPH'),...
%     get_var_value(DIRA1, set_names, var_names, 'schmid', 'schm2', 'hist2', ...
%     'PMF2', 'shapeGT2', 'fit2', 'D2', 'tMLPH'),...
%     get_var_value(DIRA2, set_names, var_names, 'schmid', 'schm2', 'hist2', ...
%     'PMF2', 'shapeGT2', 'fit2', 'D2', 'tMLPH')};
d200 = dat{1};
if any(cellfun('isempty',d200))
    return
end

% d210 = {get_var_value(DIRB0, set_names, var_names,'tau', 'wba', 'PMF2', ...
%     'shapeGT2', 'D2'),...
%     get_var_value(DIRB1, set_names, var_names, 'tau', 'wba', 'PMF2', ...
%     'shapeGT2', 'D2'),...
%     get_var_value(DIRB2, set_names, var_names, 'tau', 'wba', 'PMF2', ...
%     'shapeGT2', 'D2')};
d210 = dat{2};
if any(cellfun('isempty',d210))
    return
end

% structure data for .json format
datA0 = {};
datB0 = {};
N = numel(d200);
for n = 1:N
    
    % select results for appropriate nb. of dwell times
    dAn = d200{n};
    dBn = d210{n};
    
    % concatenate results
    datA0 = cat(1,datA0,...
        {[HA1; [dAn([1,2]),{addheadertodat(HA10,dAn{3})},...
        {addheadertodat(HA11,dAn{4})},dAn(5),{addheadertodat(HA12,dAn{6})},...
        dAn(7:8)] ]});
    datB0 = cat(1,datB0,...
        {[HB1; [dBn([1,2]),{addheadertodat(HB2,dBn{3})},dBn(4:end)] ]});
end

% transform results into a row of cell
datA0 = mat2rowcell(datA0);
datB0 = mat2rowcell(datB0);

% write file
fname = 'data_PERFONSIM1';
exportJson([dest, fname, '.json'], [H0; {datA0,datB0}], verbose);


function exportfigdat_PERFONSIM2(dat, dest, verbose)
% exportfigdat_PERFONSIM2(dat, dest, verbose)
%
% Restructure and export data for figure PERFONSIM2 of the main article.
%
% dat: {1-by-1} values of variables imported from files
% dest: destination directory
% verbose: true to show logs, false to mute

% CONSTANTS
GT_AGG_SIZES = 1:4; % D_GT=1 to 4
SAMPLE_SIZES = 1:3; % 500 to 50,000 dwell times
H0 = {'GT Db_1','GT Db_2','GT Db_3','GT Db_4'};
H1 = {'schm id','schm','hist','PMF_GT','PMF_fit','shape_GT','shape_fit',...
    'Db','computation(s)'};
H20 = {'dt','normcounts'};
H21 = {'dt','PMF'};
H22 = {'dt','fit','cumfit'};

% % show process
% if verbose
%     disp('>> collect data for figure PERFONSIM2...');
% end

% format data
num_D = length(GT_AGG_SIZES);
num_n = length(SAMPLE_SIZES);
d = cell(num_D,num_n);
for id_D = 1:num_D
    for id_n = 1:num_n
        id = sub2ind([num_n,num_D], id_n, id_D);
        d{id_D,id_n} = dat{1}{id};
    end
    if any(cellfun('isempty',d(id_D,:)))
        return
    end
end

% structure data for .json format
dat = {};
for id_n = 1:num_n
    
    % concatenate results
    dat_n = {};
    for id_D = 1:num_D
        d_Dn = d{id_D,id_n};
        dat_n = cat(2,dat_n,{[H1;[d_Dn([1,2]),{addheadertodat(H20,d_Dn{3})},...
            {addheadertodat(H21,d_Dn{4})},{addheadertodat(H22,d_Dn{5})},...
            d_Dn(6:end)]]});
    end
    dat = cat(1,dat,dat_n);
end
dat = mat2rowcell(dat);

% write file
fname = 'data_PERFONSIM2';
exportJson([dest, fname, '.json'], [H0; dat], verbose);


function exportfigdat_EBSIBS(dat, dest, verbose)
% exportfigdat_EBSIBS(dat, dest, verbose)
%
% Restructure and export data for figure EBSIBS of the main article.
%
% dat: {1-by-2} values of variables imported from files
% dest: destination directory
% verbose: true to show logs, false to mute

% CONSTANTS
H0 = {'traj','dthist_a','fit_a','BIC_a','dthist_b','fit_b','BIC_b',...
    'comp MLDPH(s)','TPM','states','pop','comp BW(s)'};
H1 = {'time','Cy3','Cy5','FRET','state'};
H2 = {'dt','hist'};
H3 = {'dt','fit','cumfit'};
H4 = {'D','schm','BIC'};

% % read data file
% if verbose
%     disp('>> collect data for figure EBSIBS...')
% end
% d1 = readData([src, FLE1], dump, verbose, 'traj');
d1 = dat{1}{1};
if isempty(d1)
    return
end
% d2 = readData([src, FLE2], dump, verbose, 'traj');
d2 = dat{1}{2};
if isempty(d2)
    return
end
% d3 = readData([src, DIRNAME], dump, verbose, 'hist1', 'fit1', 'BIC1', ...
%     'hist2', 'fit2', 'BIC2', 'tMLPH', 'TPMBW', 'states', 'pop', 'tBW');
d3 = dat{2}{1};
if isempty(d3)
    return
end

% add headers to trajectory data
traj = addheadertodat(H1, [d1{1}; d2{1}]);

% add headers to histogram data
histd1 = addheadertodat(H2, d3{1});
histd2 = addheadertodat(H2, d3{4});

% add headers to fit data
fitd1 = addheadertodat(H3, d3{2});
fitd2 = addheadertodat(H3, d3{5});

% add headers to BIC data
BIC1 = addheadertodat(H4, d3{3});
BIC2 = addheadertodat(H4, d3{6});

% structure data for .json format
dat = [{traj}, {histd1}, {fitd1}, {BIC1}, {histd2}, {fitd2}, {BIC2}, ...
    d3(7:11)];

% write file
fname = 'data_EBSIBS';
exportJson([dest, fname, '.json'], [H0; dat], verbose);


function exportfigdat_SHAPESUCCESS(dat, dest, verbose)
% exportfigdat_SHAPESUCCESS(dat, dest, verbose)
%
% Restructure and export data for figure SHAPESUCCESS of the supplementary 
% information.
%
% dat: {1-by-1} values of variables imported from files
% dest: destination directory
% verbose: true to show logs, false to mute

% defaults
GT_AGG_SIZES = 2:4; % D_GT=2 to 4
SAMPLE_SIZES = 1:3; % 500 to 50,000 dwell times
h0 = {'GT Db_2','GT Db_3','GT Db_4'};
h1 = {'schm id','schm','PMF','shape_GT','Db','spectral_GT'};
h2 = {'dt','PMF'};

% % show process
% if verbose
%     disp('>> collect data for figure SHAPESUCCESS...');
% end

% read data files
num_D = length(GT_AGG_SIZES);
num_n = length(SAMPLE_SIZES);
d = cell(num_D, num_n);
for id_D = 1:num_D
    for id_n = 1:num_n
        id = sub2ind([num_n,num_D], id_n, id_D);
        d{id_D,id_n} = dat{1}{id};
    end
    if any(cellfun('isempty',d(id_D,:)))
        return
    end
end

% structure data for .json format
dat = {};
for id_n = 1:num_n
    
    % % select results for appropriate nb. of dwell times
    % d2n = d{1,id_n};
    % d3n = d{2,id_n};
    % d4n = d{3,id_n};
    % 
    % % concatenate results
    % dat = cat(1,dat,[...
    %     {[h1;[d2n([1,2]),{addheadertodat(h2,d2n{3})},d2n(4:end)]]},...
    %     {[h1;[d3n([1,2]),{addheadertodat(h2,d3n{3})},d3n(4:end)]]},...
    %     {[h1;[d4n([1,2]),{addheadertodat(h2,d4n{3})},d4n(4:end)]]}...
    %     ]);

    % collect variables for appropriate nb. of dwell times
    dat_n = {};
    for id_D = 1:num_D
        dat_D =  d{id_D,id_n};
        dat_n = cat(2, dat_n, {[h1; [dat_D([1,2]), ...
            {addheadertodat(h2, dat_D{3})}, dat_D(4:end)]]});
    end
    
    % concatenate results
    dat = cat(1, dat, dat_n);
end

% transform results into a row of cell
dat = mat2rowcell(dat);

% write file
fname = 'data_SHAPESUCCESS';
exportJson([dest, fname, '.json'], [h0; dat], verbose);


function dat = readData(datdir,dumpdir,verbose,varargin)
% dat = readD(datdir,dumpdir,verbose,data1,data2,...)
%
% Gather analysis results through files and return them in an array.
%
% datdir: source directory OR source file for trajectory import
% dumpdir: result dump directory
% verbose: true to show logs, false to mute
% data1,data2,...: data to import from file
%   'TPMGT': {nset-by-1}[J-by-J] GT transition prob. matrix
%   'tau': [nset-by-1] GT lifetime multiplication factors
%   'wba': [nset-by-1] GT exit fraction
%   'schm[v]': {nset-by-1}[(D+2)-by-(D+2)] GT transition scheme for 
%              observed state v
%   'PMF[v]': {nset-by-1}[ndt-by-2] time vs GT PMF for observed state v
%   'shapeGT[v]': [nset-by-1] GT PMF shape code for observed state v
%   'tMLPH': [nset-by-1] ML-DPH computation time (in sec.)
%   'hist[v]': {nset-by-1}[ndt-by-2] time vs norm. histogram for observed 
%              state v
%   'fit[v]': {nset-by-1}[ndt-by-2] time vs fit PMF for observed state v
%   'BIC[v]': [nset-by-2] state degeneracy, BIC, for observed state v
%   'D[v]': [nset-by-1] inferred state degeneracy for observed state v
%   'traj': [L-by-5] time, don and acc intensity, FRET and state 
%           trajectories
%   'TPMBW': {nset-by-1}[J-by-J] inferred transition prob. matrix
%   'tBW': [nset-by-1] BW computation time (in sec.)
% dat: {1-by-nDat} imported data

% defaults
DTMAX = 1E9;
PMIN = 1E-9;
NDT = 1000;
SPLT0 = 0.1; % reference trajectory sampling time in simulated data
EXT_DATA = '_data.mat';
EXT_SIMPRM = '_simprm.mat';
EXT_SIMRES = '_simres.mat';
EXT_DPH = '_mldphres.mat';
EXT_BW = '_bwres.mat';

% initialize output
dat = [];

% collect histogram binning for experimental data
dphprm = PHtest_adjustparam(PHtest_getdefaultparam(),'EBS-IBS');
bin_exp = dphprm.bin;
excl_exp = dphprm.excl;

% get dump folder specific to analysis
dumpdir = PHtest_getdumpflddir(dumpdir,datdir);
if ~isempty(dumpdir) 
    if dumpdir(end)~=filesep
        dumpdir = [dumpdir,filesep];
    end
end

% check if it's only about collecting trajectories from a single file
if any(contains(varargin,{'traj'}))

    % only import trajectories from a single file
    if ~exist(datdir,'file')
        if verbose
            disp(['trajectory file not found: ',datdir]);
        end
        return
    end

    [~,~,fext] = fileparts(datdir);
    switch fext
        case '.txt'
            fretraj = importdata(datdir,'\t',3);
            traj = fretraj.data(:,[1,3,4,7,8]);
        case '.traces'
            fretraj = importdata(datdir,'\t',3);
            traj = fretraj.data(:,[1,5:8]);
        case '.mat'
            simdat = load(datdir,'res');
            dt = simdat.res.dt_gt;
            N = numel(dt);
            traj = cell(N,1);
            for n = 1:N
                traj{n} = getDiscrFromDt(dt{n}(:,[1,2]),SPLT0);
            end
    end
    dat = {{traj}};
    return
end

% check souce directory
if datdir(end)~=filesep
    datdir = [datdir,filesep];
end
if ~exist(datdir,'dir')
    if verbose
        disp(['folder ',datdir,' not found.']);
    end
    return
end

% check whether source directory contains simulated or experimental
% data based on the presence of simulation parameter files
flist_prm = dir([datdir,'*',EXT_SIMPRM]);
if isempty(flist_prm)
    [~,subdir,~] = fileparts(datdir(1:end-1));
    flist_prm = dir([datdir,'*',EXT_DATA]);
    if isempty(flist_prm) && verbose
        disp(['Neither simulation parameters or experimental data were ',...
            'found in sub-directory ',subdir])
        return
    else
        issim = false; % experiment: fset contains the exp. data file
    end
else
    issim = true; % simulation: fset lists parameter files
    F1 = size(flist_prm,1);
end

% list all interresting files
fprm = {};
fdat = {};
fres = {};
fbw = {};
bin = [];
excl = [];
if verbose
    disp('List necessary files ...');
end
if issim
    % collect for analysis of simulated data
    for f1 = 1:F1

        % list simulated data sets
        fprm_n = [datdir,flist_prm(f1,1).name];
        setname = flist_prm(f1).name(1:end-length(EXT_SIMPRM));
        if ~exist(fprm_n,'file')
            if verbose
                fprintf(' no parmeter file for set %s.\n',setname);
            end
        end

        % list simulated data files
        flist_dat = dir([datdir,'*',setname,'*',EXT_SIMRES]);
        F2 = size(flist_dat,1);
        src_res = [datdir,dumpdir,setname,filesep];
        for f2 = 1:F2
            fdat_n = [flist_dat(f2).folder,filesep,flist_dat(f2).name];
            if ~exist(fdat_n,'file')
                fdat_n = {''};
            else
                dat0 = struct2table(whos('-file',fdat_n));
                if isempty(dat0.name) || ~any(contains(dat0.name,'res'))
                    fdat_n = {''};
                end
            end
            
            datname = flist_dat(f2).name(1:end-length(EXT_SIMRES));
            fres_n = [src_res,datname,EXT_DPH];
            if ~exist(fres_n,'file')
                fres_n = {''};
            end
                
            fbw_n = [src_res,datname,EXT_BW];
            if ~exist(fbw_n,'file')
                fbw_n = {''};
            end

            fprm = cat(1,fprm,fprm_n);
            fdat = cat(1,fdat,fdat_n);
            fres = cat(1,fres,fres_n);
            fbw = cat(1,fbw,fbw_n);
            bin = cat(1,bin,1);
            excl = cat(1,excl,false);
        end
    end

else
    % collect for analysis of experimental data
    [~,srcdir,~] = fileparts(flist_prm.folder);
    if verbose
        disp(['process file ',[srcdir,filesep,flist_prm.name],'...']);
    end

    fdat_n = [flist_prm.folder,filesep,flist_prm.name];
    if ~exist(fdat_n,'file')
        fdat_n = {''};
    end
    
    % list ML-PH result files
    flist_res = dir([datdir,dumpdir,'*',EXT_DPH]);
    if size(flist_res,1)>0
        fres_n = [flist_res.folder,filesep,flist_res.name];
        if ~exist(fres_n,'file')
            fres_n = {''};
        end
        
        datname = flist_res.name(1:end-length(EXT_DPH));
        fbw_n = [flist_res.folder,filesep,datname,EXT_BW];
        if ~exist(fbw_n,'file')
            fbw_n = {''};
        end
    else
        fres_n = {''};
        fbw_n = {''};
    end

    fprm = cat(1,fprm,{''});
    fdat = cat(1,fdat,fdat_n);
    fres = cat(1,fres,fres_n);
    fbw = cat(1,fbw,fbw_n);
    bin = cat(1,bin,bin_exp);
    excl = cat(1,excl,excl_exp);
end

% determine dwell time axis for GT and fit plot
xgt = unique(round(logspace(0,log10(DTMAX),NDT)));

% collect data from files
N = size(fdat,1);
ndat = numel(varargin);
dat = cell(1,ndat);
for n = 1:N
    if isempty(fdat{n})
        continue
    end
    if verbose
        dprm = fileparts(fprm{n});
        [~,srcdir,~] = fileparts(dprm);
        [~,namedat,extdat] = fileparts(fdat{n});
        disp([sprintf('process file %*i/%i: ',nbdigit(N),n,N),...
            [srcdir,filesep,namedat,extdat],'...']);
    end

    filedata_n = init_filedata(1);

    if ~isempty(fprm{n})
        % load simulation parameters
        prm = load(fprm{n});
        filedata_n(1,"prm") = {{prm}};
            
        % preliminary calculations
        if any(contains(varargin,{'schm','PMF','shapeGT','spectralGT'}))
            [T,t,ip,schm,~] = calc_absmm_from_hmm(prm.ip,prm.tp,prm.val);
            [fprm_src,~,~] = fileparts(fprm{n});
            setname = fprm{n}((length(fprm_src)+2):(end-length(EXT_SIMPRM)));
            schmid = split(setname,'_');
            schmid = str2double(schmid{end});
            filedata_n(1,["schm","schmid"]) = {{schm},{schmid}};
        end
        if any(contains(varargin,'tau'))
            tau0 = (1./sum(prm.tp,2))';
            taub0 = tau0(prm.val==prm.val(end));
            fact = sort(taub0);
            fact = round(fact(2)/fact(1));
            filedata_n(1,"fact") = {{fact}};
        end
        if any(contains(varargin,'wba'))
            Da0 = sum(prm.val==prm.val(1));
            wba0 = prm.tp((Da0+1):end,1)./sum(prm.tp((Da0+1):end,:),2);
            wba0 = wba0(1);
            filedata_n(1,"wba0") = {{wba0}};
        end
        if any(contains(varargin,'shapeGT'))
            [~,~,~,~,code] = effectivedegeneracy_from_dlogP_dlogt(T,ip,xgt,...
                PMIN);
            filedata_n(1,"codeGT") = {{code}};
        end
        if any(contains(varargin,'spectralGT'))
            V = length(T);
            spectral = cell(1,V);
            for v = 1:V
                [a, eigval] = calcexpweight(T{v}, ip{v}, 0, 1);
                spectral{v} = [a(:)'; eigval(:)'];
            end
            filedata_n(1,"spectralGT") = {{spectral}};
        end
    end
    if any(contains(varargin,{'hist','PMF','fit'}))
        if issim
            dat0 = load(fdat{n},'res');
            dat0 = dat0.res;
        else
            [prm,dat0] = PHtest_importexpdata(fdat{n},'');
        end
        [dtdat,Pdat,edgdat] = ...
            builddthist(dat0.dt_obs,prm.val,prm.rate,bin(n),excl(n));
        V = numel(unique(prm.val));
        xdat = cell(1,V);
        for v = 1:V
            xdat{v} = dtdat{v}(Pdat{v}>0); % unique(dt_obs)
        end
        filedata_n(1,["edgdat","xdat","Pdat"]) = {{edgdat},{xdat},{Pdat}};

        if any(contains(varargin,{'PMF','fit'})) && ~isempty(fprm{n})
            Pgt = cell(1,V);
            for v = 1:numel(T)
                Pgt{v} = calc_DPH_PMF(T{v},t{v},ip{v},xgt);
                % xgt1 = edgdat{v}(1:end-1)-1; % [0,1,2,...,(dt_obs_max-1)]
                % xgt1 = xgt1(Pdat{v}>0);
                % xgt2 = edgdat{v}(2:end)-1; % [1,2,...,dt_obs_max]
                % xgt2 = xgt2(Pdat{v}>0);
                % CDF1 = calcDPHCDF(T{v},ip{v},xgt1,0);
                % CDF2 = calcDPHCDF(T{v},ip{v},xgt2,0);
                % Pgt{v} = CDF2-CDF1;
            end
            filedata_n(1,["xgt","Pgt"]) = {{xgt},{Pgt}};
        end
    end
    if ~isempty(fres{n})
        res = load(fres{n},'-mat'); 
        if isfield(res,'phprm')
            analysismethod = ...
                PHtest_getmethodfromcalcmode(res.phprm.calcmode);
        else
            analysismethod = 'mlph';
        end
        res = res.dphres;
        filedata_n(1,["analysismethod","res"]) = {analysismethod,{res}};
    end
    if ~isempty(fbw{n})
        bw = load(fbw{n},'bwres_w'); 
        bw = bw.bwres_w;
        filedata_n(1,"bw") = {{bw}};
    end

    % pick and format requested data
    for d = 1:ndat
        dat{d} = cat(1,dat{d},...
            pick_data_from_res(varargin{d},issim,filedata_n,PMIN,verbose));
    end
end


function filedata = init_filedata(N)
filedata.prm = {};
filedata.fact = [];
filedata.wba0 = [];
filedata.schm = {};
filedata.schmid = [];
filedata.xdat = {};
filedata.edgdat = {};
filedata.xgt = {};
filedata.Pgt = {};
filedata.Pdat = {};
filedata.codeGT = {};
filedata.spectralGT = {};
filedata.res = {};
filedata.bw = {};
filedata.analysismethod = {};

% trick: add dummy line because repmat returns empty table for N=1
filedata = struct2table(repmat(filedata,N+1,1)); 
filedata = filedata(1:N,:); % remove dummy line


function dat = pick_data_from_res(name,issim,fdat,pmin,verbose)
dat = [];
switch name
    case {'TPMGT','tau','wba','schm1','schm2','schmid','PMF1','PMF2',...
            'shapeGT1','shapeGT2','spectralGT1','spectralGT2'}

        if ~issim 
            if verbose
                disp(['PHtest_collectResults>checksim: unknown data ',...
                    '\"',name,'\" for experimental set.']);
            end
            return
        end

        switch name
            case 'TPMGT'
                dat = fdat.prm{1}.tp;
            case 'tau'
                dat = fdat.fact;
            case 'wba'
                dat = fdat.wba0;
            case {'schm1','schm2'}
                v = str2double(name(end));
                dat = fdat.schm{1}(v);
            case 'schmid'
                dat = fdat.schmid{1};
            case {'PMF1','PMF2'}
                v = str2double(name(end));
                dat = {[fdat.xgt{1}',fdat.Pgt{1}{v}']};
            case {'shapeGT1','shapeGT2'}
                v = str2double(name(end));
                dat = fdat.codeGT{1}(v);
            case {'spectralGT1','spectralGT2'}
                v = str2double(name(end));
                dat = fdat.spectralGT{1}(v);
        end

    case {'hist1','hist2'}
        v = str2double(name(end));
        dat = {[fdat.xdat{1}{v}',fdat.Pdat{1}{v}(fdat.Pdat{1}{v}>0)']};
    
    case 'tMLPH'
        if ~isempty(fdat.res{1})
            switch fdat.analysismethod{1}
                case 'mlph'
                    dat = fdat.res{1}{3}.t_dphtest;
                case 'emexp'
                    dat = fdat.res{1}{3}.t_emexp;
                case 'iemm'
                    dat = fdat.res{1}{3}.t_iemm;
                case 'iamm'
                    dat = fdat.res{1}{3}.t_iamm;
            end
        else
            dat = NaN;
        end

    case {'fit1','fit2'}
        v = str2double(name(end));
        if ~isempty(fdat.res{1})
            switch fdat.analysismethod{1}
                case {'mlph','iamm'}
                    Tfit = fdat.res{1}{3}.tp_fit{v}(:,1:end-1);
                    try
                        tfit = fdat.res{1}{3}.tp_fit{v}(:,end);
                    catch err
                        throw(err);
                    end
                    ipfit = fdat.res{1}{3}.pi_fit{v};     
                    if issim
                        xfit = fdat.xgt{1};
                        Pfit = calc_DPH_PMF(Tfit,tfit,ipfit,xfit);
                        cumPfit = cumsum(Pfit);
                    else
                        xfit = fdat.xdat{1}{v};
                        xfit1 = fdat.edgdat{1}{v}(1:end-1)-1;
                        xfit1 = xfit1(fdat.Pdat{1}{v}>0);
                        xfit2 = fdat.edgdat{1}{v}(2:end)-1;
                        xfit2 = xfit2(fdat.Pdat{1}{v}>0);
                        CDF1 = calcDPHCDF(Tfit,ipfit,xfit1,eps("double"));
                        CDF2 = calcDPHCDF(Tfit,ipfit,xfit2,eps("double"));
                        Pfit = CDF2-CDF1;
                        cumPfit = CDF1;
                    end
                    dat = {[xfit',Pfit',cumPfit']};

                case {'emexp','iemm'}
                    if issim
                        xfit = fdat.xgt{1};
                    else
                        xfit = fdat.xdat{1}{v};
                    end
                    taufit = fdat.res{1}{3}.tau{v};
                    afit = fdat.res{1}{3}.a{v};
                    Pfit = zeros(size(xfit));
                    for d = 1:numel(afit)
                        Pfit = Pfit + afit(d)*exppdf(xfit,taufit(d));
                    end
                    % Pfit = sum(repmat(a(:),1,numel(xfit)).*...
                    %     geom_density(taufit,xfit),1);
                    cumPfit = cumsum(Pfit/sum(Pfit));
                    % dat = {[xfit',Pfit'/sum(Pfit),cumPfit']};
                    dat = {[xfit',Pfit',cumPfit']};
            end
        else
            dat = {[NaN,NaN,NaN]};
        end

    case {'shapefit1','shapefit2'}
        v = str2double(name(end));
        if ~isempty(fdat.res{1})  
            if issim
                xfit = fdat.xgt{1};
            else
                xfit = fdat.xdat{1}{v};
            end
            switch fdat.analysismethod{1}
                case {'mlph','iamm'}
                    Tfit = fdat.res{1}{3}.tp_fit{v}(:,1:end-1);
                    ipfit = fdat.res{1}{3}.pi_fit{v}; 
                    [~,~,~,~,code] = effectivedegeneracy_from_dlogP_dlogt(...
                        Tfit,ipfit,xfit,pmin);  

                case {'emexp','iemm'}
                    Dfit = length(fdat.res{1}{3}.tau{v});
                    code = 2*ones(1,Dfit-1);
                    if isempty(code)
                        code = {[]};
                    end
            end
            dat = {code};
        else
            dat = {NaN};
        end

    case {'BIC1','BIC2'}
        v = str2double(name(end));
        if ~isempty(fdat.res{1})
            dat = {getBICres(fdat.res{1}{5}{v})};
        else
            dat = {cell(1,3)};
        end

    case {'D1','D2'}
        v = str2double(name(end));
        if ~isempty(fdat.res{1})
            switch fdat.analysismethod{1}
                case {'mlph','iamm'}
                    dat = size(fdat.res{1}{3}.schm{v},1)-2;

                case {'emexp','iemm'}
                    dat = fdat.res{1}{3}.D(v);
            end
        else
            dat = NaN;
        end

    case 'states'
        if ~isempty(fdat.res{1})
            dat = fdat.res{1}(4);
        else
            dat = {};
        end

    case 'TPMBW'
        if ~isempty(fdat.bw{1})
            dat = {cat(3,fdat.bw{1}{1},fdat.bw{1}{2})};
        else
            dat = {};
        end

    case 'pop'
        if ~isempty(fdat.bw{1})
            states = unique(fdat.bw{1}{4}.dt(:,3))';
            pop = zeros(1,max(states));
            for j = states
                pop(j) = sum(fdat.bw{1}{4}.dt(fdat.bw{1}{4}.dt(:,3)==j,1));
            end
            pop = pop/sum(pop);
            dat = {pop};
        else
            dat = {};
        end

    case 'tBW'
        if ~isempty(fdat.bw{1})
            dat = fdat.bw{1}{5};
        else
            dat = NaN;
        end

    otherwise
        if verbose
            disp(['PHtest_collectResults>readData: unknown '...
                'data \"',name,'\".'])
        end
end


function rowcell = mat2rowcell(mat)
% rowcell = mat2rowcell(mat)
%
% Converts a matrix into a row cell array, each matrix column being stored
% in one cell
%
% mat: [R-by-C] matrix
% rowcell: {1-by-C} cell array

[~,C] = size(mat);
rowcell = cell(1,C);
for c = 1:C
    rowcell{1,c} = mat(:,c);
end


function BICres = getBICres(mdl)
S = size(mdl,1);
BICres = cell(1,3);
BICres{1} = NaN(S,1);
BICres{2} = cell(S,1);
BICres{3} = NaN(S,1);
for s = 1:S
    BICres{1}(s) = size(mdl(s).schm,1)-2;
    BICres{2}{s} = mdl(s).schm;
    BICres{3}(s) = mdl(s).BIC;
end


function dat = addheadertodat(hd,dat0)
% dat = addheadertodat(hd,dat0)
%
% Add headers to each each cell in dat0.
%
% dat0: {S-by-1}[ns-by-C] data
% hd: {1-by-C} headers
% dat: {2-by-C} data with headers

dat = {};
for s = 1:numel(dat0)
    dat = cat(1,dat,{[hd;mat2rowcell(dat0{s})]});
end


function traj = getFirstTrajStartingWithJ(dat,j)
% traj = getFirstTrajStartingWithJ(dat,j)
%
% Rreturns state trajectories of the first molecule in the input initiated
% with states j.
%
% dat: {1-by-N}[Ln-by-2] molecule time and state trajectories
% j: [1-by-J] starting states
% traj: [L-by-2] selected molecule time and state trajectories

N = numel(dat);
n = 1;
while n<N && ~any(dat{n}(1)==j)
    n = n+1;
end
traj = dat{n};


function [traj,n] = getFirstTrajContainingJ(dat,j,minfrac)
% [traj,n] = getFirstTrajContainingJ(dat,j,minfrac)
%
% Returns state trajectories of the first molecule in the input containing 
% all states j.
%
% dat: {1-by-N}[Ln-by-2] molecule time and state trajectories
% j: [1-by-J] states
% traj: [L-by-2] selected molecule time and state trajectories
% n: selected molecule index

N = numel(dat);
n = 1;
while n<N && ~all(sum(repmat(dat{n},1,size(j,2))==repmat(j,size(dat{n},1),...
        1),1)/size(dat{n},1)>minfrac)
    n = n+1;
end
traj = dat{n};


function issim = checksim(issim,dataname,verbose)
if ~issim && verbose
    disp(['PHtest_collectResults>checksim: unknown data ',...
        '\"',dataname,'\" for experimental set.']);
end


function exportJson(jfile,dat,verbose)
% exportJson(jfile,dat,verbose)
%
% Export data into a .json file using JSON formatting
%
% jfile: destination file address
% dat: {2-by-C} header strings and numerical data columns
% verbose: true to show logs, false to mute

[pname,~,~] = fileparts(jfile);
if ~exist(pname,'dir')
    mkdir(pname);
end

f = fopen(jfile,'w');
writeJson(f,dat,0);
fclose(f);
if verbose
    disp(['File ',jfile,' was successfully exported!']);
end


function writeJson(f,dat,ntab)
% writeJson(f,dat,ntab)
%
% Write header and data into an ASCII file using JSON formatting
%
% f: file identifier
% dat: {2-by-C} header strings (1st row) and associated numerical data (2nd 
%  row) or [R-by-C] numerical data.
% ntab: data indentation

[R,C] = size(dat);
if R==0 || C==0
    fprintf(f,'[]');
    return
end

tabs = repmat('\t',[1,ntab]);
ishd = iscell(dat) && R==2 && ...
    (all(cellfun(@ischar,dat(1,:)) | cellfun(@isempty,dat(1,:))));

if C~=1 
    if ishd
        fprintf(f,'{');
    else
        fprintf(f,'[');
    end
end
for c = 1:C
    fprintf(f,['\n',tabs,'\t']);
    if ishd
        if ~isempty(dat{1,c}) && ischar(dat{1,c}) % write header
            fprintf(f,['\"',dat{1,c},'\": ']);
        end
        if iscell(dat{2,c})
            writeJson(f,dat{2,c},ntab+1);
        else
            [Rc,Cc] = size(dat{2,c});
            if Rc~=1
                fprintf(f,'[');
            end
            for r = 1:Rc
                if Cc>0
                    fmtdat = '%d';
                    if Cc>1
                        fprintf(f,'[');
                        fmtdat = cat(2,fmtdat,repmat(',%d',[1,Cc-1]));
                    end
                    fprintf(f,fmtdat,dat{2,c}(r,:));
                    if Cc>1
                        fprintf(f,']');
                    end
                end
                if r~=Rc
                    fprintf(f,',');
                end
            end
            if Rc~=1
                fprintf(f,']');
            end
        end

    else
        if R~=1
            fprintf(f,'[');
        end
        for r = 1:R
            if iscell(dat{r,c})
                writeJson(f,dat{r,c},ntab+1);
            else
                [Rc,Cc] = size(dat{r,c});
                if Rc~=1
                    fprintf(f,'[');
                end
                for rc = 1:Rc
                    if Cc>0
                        fmtdat = '%d';
                        if Cc>1
                            fprintf(f,'[');
                            fmtdat = cat(2,fmtdat,repmat(',%d',[1,Cc-1]));
                        end
                        fprintf(f,fmtdat,dat{r,c}(rc,:));
                        if Cc>1
                            fprintf(f,']');
                        end
                    end
                    if rc~=Rc
                        fprintf(f,',');
                    end
                end
                if Rc~=1
                    fprintf(f,']');
                end
            end
            if r~=R
                fprintf(f,',');
            end
        end
        if R~=1
            fprintf(f,']');
        end
    end

    if c~=C
        fprintf(f,',');
    end
end
if C~=1
    if ishd
        fprintf(f,['\n',tabs,'}']);
    else
        fprintf(f,['\n',tabs,']']);
    end
end


