function DPHtest_collectResults(rootdir,destdir,varargin)
% DPHtest_collectResults(rootdir,destdir)
% DPHtest_collectResults(rootdir,destdir,subdir)
%
% Collect results for plotting figures and epxort them to .json file format.
%
% rootdir: source directory
% destdir: destination directory
% subdir: (opt) {1-by-D} specific subdirectories to export

% set MATLAB search path
codePath = fileparts(mfilename('fullpath'));
addpath(genpath(codePath));

% check existence of source directory
if ~exist(rootdir,'dir')
	disp('Source directory not found.')
	return
end
if rootdir(end)~=filesep
	rootdir = [rootdir,filesep];
end

% create destination directory if not existing
if ~exist(destdir,'dir')
    disp('create destination directory...')
    mkdir(destdir);
end
if destdir(end)~=filesep
	destdir = [destdir,filesep];
end

% collect subdirectories
if ~isempty(varargin)
    subfig = varargin{1};
    if ~iscell(subfig)
        disp(['DPHtest_collectResults: ERROR argument 3 must be a ',...
            'cell array']);
        return
    end
else
    subfig = {'figure2A','figure2B','figure2BC','figure2C','figure2D',...
        'figure3A','figure3B','figure4B','figure4C','figure4D','figure4E'};
end
strsubfig = '';
for sd = 1:numel(subfig)
    strsubfig = cat(2,strsubfig,subfig{sd},' ');
end
disp(['files will be exported for: ',strsubfig(1:end-1),'.']);

% figure 2A: D and model accuracy = f(D_GT) for models I, II and III
if any(cellstrcmp(subfig,'figure2A'))
    disp('process subfigure 2A...')
    exportfigdat_2A(rootdir,destdir,false);
    exportfigdat_2A(rootdir,destdir,true);
end

% figure 2B: D and model accuracy = f(tau_b2) for model I and II
if any(cellstrcmp(subfig,'figure2B'))
    disp('process subfigure 2B...')
    exportfigdat_2B(rootdir,destdir,false);
    exportfigdat_2B(rootdir,destdir,true);
end

% figure 2BC: D = f(tau_b2,Nab/Nbb) for model II and D=2
if any(cellstrcmp(subfig,'figure2BC'))
    disp('process subfigure 2BC...')
    exportfigdat_2BC(rootdir,destdir,false);
    exportfigdat_2BC(rootdir,destdir,true);
end

% figure 2C: D and model accuracy = f(Nab/Nbb) for model II and three DPH fit curves
if any(cellstrcmp(subfig,'figure2C'))
    disp('process subfigure 2C...')
    exportfigdat_2C(rootdir,destdir,false);
    exportfigdat_2C(rootdir,destdir,true);
end

% figure 2D: D and model accuracy = f(state connexions) 
if any(cellstrcmp(subfig,'figure2D'))
    disp('process subfigure 2D...')
    exportfigdat_2D(rootdir,destdir,false);
    exportfigdat_2D(rootdir,destdir,true);
end

% figure 3A: GT and inferred TPMs for quenched disorder, ML-DPH BICs and TPs, and computation time
if any(cellstrcmp(subfig,'figure3A'))
    disp('process subfigure 3A...')
    exportfigdat_3A(rootdir,destdir,false)
    exportfigdat_3A(rootdir,destdir,true)
end

% figure 3B: GT and inferred TPM for dynamic disorder, ML-DPH BICs and TPs and computation time
if any(cellstrcmp(subfig,'figure3B'))
    disp('process subfigure 3B...')
    exportfigdat_3B(rootdir,destdir,false);
    exportfigdat_3B(rootdir,destdir,true);
end

% figure 4B: EBS-IBS intensity- and FRET-time traces
if any(cellstrcmp(subfig,'figure4B'))
    disp('process subfigure 4B...')
    datdir = 'EBS-IBS/trajectories/IBS-20mM-Mg-3_all185_mol113_post(2).traces';
    hdrs_4A = {'time','Cy3','Cy5','FRET','state'};
    dat = readData([rootdir,datdir],'traj');
    dat_4A = mat2rowcell(dat);
    exportJson([destdir,'figure4',filesep,'data_figure4B.json'],...
        [hdrs_4A;dat_4A]);
end

% figure 4CDE: EBS-IBS ML-DPH BICs, dwell time histograms best fit curves, 
% inferred TPM and computation time
if any(cellstrcmp(subfig,'figure4CDE'))
    disp('process subfigure 4CDE...')
    datdir = 'EBS-IBS';
    
    hdrs_4C11 = {'D','scheme','BIC_a','BIC_b','std(BIC_a)','std(BIC_b)'};
    hdrs_4C12 = {'STPM_a','std(STPM_a)','STPM_b','std(STPM_b)'};
    hdrs_4C13 = {'t','std(t)','N','sampling (s)'};
    hdrs_4C145 = {'dt','mean count','std count','mean fit','std fit'};
    hdrs_4C22 = {'t','std(t)'};
    hdrs_4C21 = {'states','TPM_wo','-err_wo','+err_wo','std(TPM_wo)',...
        'std(-err_wo)','std(+err_wo)'};
    hdrs_4C1 = {'BIC','STPM','t_comp','dthista','dthistb'};
    hdrs_4C2 = {'TPM','t_comp'};
    hdrs_4C = {'DPH','BW'};
    
    disp('>> collect data for DPH fit model...')
    dat = readData([rootdir,datdir],'model');
    dat_4C11 = [hdrs_4C11; mat2rowcell(dat{2})];
    dat_4C12 = [hdrs_4C12; dat{3}];
    dat_4C13 = [hdrs_4C13; dat{4}];
    dat_4C14 = [hdrs_4C145; mat2rowcell(dat{5}{1})];
    dat_4C15 = [hdrs_4C145; mat2rowcell(dat{5}{2})];
    dat_4C1 = [hdrs_4C1;{dat_4C11,dat_4C12,dat_4C13,dat_4C14,dat_4C15}];
    dat_4C21 = [hdrs_4C21;dat{1}];
    dat_4C22 = [hdrs_4C22;dat{6}];
    dat_4C2 = [hdrs_4C2;{dat_4C21,dat_4C22}];
    dat_4C = {dat_4C1,dat_4C2};
    exportJson([destdir,'figure4',filesep,'data_figure4CDE.json'],...
        [hdrs_4C;dat_4C]);
end

disp('Process completed!');


function exportfigdat_2A(rootdir,destdir,sumexp)

% defaults
datdir_I = 'dataset1';
datdir_II = 'dataset2';
datdir_III = 'dataset3';
hdrs = {'GT Db_I','Db_I','std(Db_I)','acc_I','std(acc_I)',...
    'GT Db_II','Db_II','std(Db_II)','acc_II','std(acc_II)',...
    'GT Db_III','Db_III','std(Db_III)','acc_III','std(acc_III)',...
    'dthistI','dthistII','dthistIII'};
hdrs2 = {'dt','mean count','std count','mean fit','std fit','GT'};

% show process
if sumexp
    disp('>> collect data for sum exp. fit model...')
else
    disp('>> collect data for DPH fit model...');
end

% read data files
dat_I = readData([rootdir,datdir_I],'none',sumexp);
dat_II = readData([rootdir,datdir_II],'none',sumexp);
dat_III = readData([rootdir,datdir_III],'none',sumexp);

% structure data for .json format
dat1 = [dat_I{1},dat_II{1},dat_III{1}];
dat1 = mat2rowcell(dat1);
dat2_I = {};
dat2_II = {};
dat2_III = {};
for n = 1:numel(dat_I{2})
    dat2_I = cat(1,dat2_I,{[hdrs2;mat2rowcell(dat_I{2}{n})]});
end
for n = 1:numel(dat_II{2})
    dat2_II = cat(1,dat2_II,{[hdrs2;mat2rowcell(dat_II{2}{n})]});
end
for n = 1:numel(dat_III{2})
    dat2_III = cat(1,dat2_III,{[hdrs2;mat2rowcell(dat_III{2}{n})]});
end
dat = [dat1,{dat2_I},{dat2_II},{dat2_III}];

% write file
fname = 'data_figure2A';
if sumexp
    fname = [fname,'_sumexp'];
end
exportJson([destdir,'figure2',filesep,fname,'.json'],[hdrs;dat]);


function exportfigdat_2B(rootdir,destdir,sumexp)

% defaults
datdir_I = 'dataset5';
datdir_II = 'dataset6';
hdrs = {'GT Db_I','GT tau_b1_I','GT tau_b2_I','Db_I','std(Db_I)',...
    'acc_I','std(acc_I)','GT Db_II','GT tau_b1_II','GT tau_b2_II',...
    'Db_II','std(Db_II)','acc_II','std(acc_II)','dthistI','dthistII'};
hdrs2 = {'dt','mean count','std count','mean fit','std fit','GT'};

% show process
if sumexp
    disp('>> collect data for sum exp. fit model...')
else
    disp('>> collect data for DPH fit model...');
end

% read data files
dat_I = readData([rootdir,datdir_I],'tau',sumexp);
dat_II = readData([rootdir,datdir_II],'tau',sumexp);

% structure data for .json format
dat1 = [dat_I{1},dat_II{1}];
dat1 = mat2rowcell(dat1);
dat2_I = {};
for n = 1:numel(dat_I{2})
    dat2_I = cat(1,dat2_I,{[hdrs2;mat2rowcell(dat_I{2}{n})]});
end
dat2_II = {};
for n = 1:numel(dat_II{2})
    dat2_II = cat(1,dat2_II,{[hdrs2;mat2rowcell(dat_II{2}{n})]});
end
dat = [dat1,{dat2_I},{dat2_II}];

% write file
fname = 'data_figure2B';
if sumexp
    fname = [fname,'_sumexp'];
end
exportJson([destdir,'figure2',filesep,fname,'.json'],[hdrs;dat]);


function exportfigdat_2BC(rootdir,destdir,sumexp)

% defaults
datdir = 'dataset10';
hdrs = {'GT Db','GT fact','GT tau_b1','GT tau_b2','GT Nab/Nbb','Db',...
    'std(Db)','dthist'};
hdrs2 = {'dt','mean count','std count','mean fit','std fit','GT'};

% show process
if sumexp
    disp('>> collect data for sum exp. fit model...')
else
    disp('>> collect data for DPH fit model...');
end

% read data files
dat = readData([rootdir,datdir],'tauwba',sumexp);

% structure data for .json format
dat1 = mat2rowcell(dat{1});
dat2 = {};
for n = 1:numel(dat{2})
    dat2 = cat(1,dat2,{[hdrs2;mat2rowcell(dat{2}{n})]});
end
dat = [dat1,{dat2}];

% write file
fname = 'data_figure2BC';
if sumexp
    fname = [fname,'_sumexp'];
end
exportJson([destdir,'figure2',filesep,fname,'.json'],[hdrs;dat]);


function exportfigdat_2C(rootdir,destdir,sumexp)

% defaults
datdir = 'dataset4';
hdrs = {'GT Db_II','Nab/Nbb','Db_II','std(Db_II)','acc_II','std(acc_II)',...
    'dthist'};
hdrs2 = {'dt','mean count','std count','mean fit','std fit','GT'};

% show process
if sumexp
    disp('>> collect data for sum exp. fit model...')
else
    disp('>> collect data for DPH fit model...');
end

% read data files
dat = readData([rootdir,datdir],'wba',sumexp);

% structure data for .json format
dat1 = mat2rowcell(dat{1});
dat2 = {};
for n = 1:numel(dat{2})
    dat2 = cat(1,dat2,{[hdrs2;mat2rowcell(dat{2}{n})]});
end
dat = [dat1,{dat2}];

% write file
fname = 'data_figure2C';
if sumexp
    fname = [fname,'_sumexp'];
end
exportJson([destdir,'figure2',filesep,fname,'.json'],[hdrs;dat]);


function exportfigdat_2D(rootdir,destdir,sumexp)

% defaults
datdir = 'dataset9';
hdrs = {'GT Db','GT schm','Db','std(Db)','acc','std(acc)','dthist'};
hdrs2 = {'dt','mean count','std count','mean fit','std fit','GT'};

% show process
if sumexp
    disp('>> collect data for sum exp. fit model...')
else
    disp('>> collect data for DPH fit model...');
end

% read data files
dat = readData([rootdir,datdir],'schm',sumexp);

% structure data for .json format
dat1 = mat2rowcell(dat{1});
dat2 = {};
for n = 1:numel(dat{3})
    dat2 = cat(1,dat2,{[hdrs2;mat2rowcell(dat{3}{n})]});
end
dat = [dat1(1),dat(2),dat1(2:end),{dat2}];

% write file
fname = 'data_figure2D';
if sumexp
    fname = [fname,'_sumexp'];
end
exportJson([destdir,'figure2',filesep,fname,'.json'],[hdrs;dat]);


function exportfigdat_3A(rootdir,destdir,sumexp)

% defaults
datdir = 'dataset7';
hdrs = {'DPH','BW'};
hdrs1 = {'BIC','STPM','t_comp','dthista','dthistb'};
hdrs2 = {'TPM','t_comp'};
hdrs11 = {'D','scheme','BIC_a','BIC_b','std(BIC_a)','std(BIC_b)'};
hdrs12 = {'STPM_a','std(STPM_a)','STPM_b','std(STPM_b)'};
hdrs13 = {'t','std(t)','N','sampling (s)'};
hdrs14 = {'dt','mean count','std count','mean fit','std fit','GT'};
hdrs22 = {'t','std(t)'};
hdrs21 = {'GT states','GT TPM','states','TPM_w','-err_w','+err_w','TPM_wo',...
    '-err_wo','+err_wo','std(TPM_w)','std(-err_w)','std(+err_w)',...
    'std(TPM_wo)','std(-err_wo)','std(+err_wo)'};

% show process
if sumexp
    disp('>> collect data for sum exp. fit model...')
else
    disp('>> collect data for DPH fit model...');
end

% read data files
dat = readData([rootdir,datdir],'model',sumexp);

% structure data for .json format
dat11 = [hdrs11; mat2rowcell(dat{2})];
dat12 = [hdrs12; dat{3}];
dat13 = [hdrs13; dat{4}];
dat14 = [hdrs14; mat2rowcell(dat{5}{1})];
dat15 = [hdrs14; mat2rowcell(dat{5}{2})];
dat1 = [hdrs1; {dat11,dat12,dat13,dat14,dat15}];
dat21 = [hdrs21; dat{1}];
dat22 = [hdrs22; dat{6}];
dat2 = [hdrs2; {dat21,dat22}];
dat = {dat1,dat2};

% write file
fname = 'data_figure3A';
if sumexp
    fname = [fname,'_sumexp'];
end
exportJson([destdir,'figure3',filesep,fname,'.json'],[hdrs;dat]);


function exportfigdat_3B(rootdir,destdir,sumexp)

% defaults
datdir = 'dataset8';
hdrs11 = {'D','scheme','BIC_a','BIC_b','std(BIC_a)','std(BIC_b)'};
hdrs12 = {'STPM_a','std(STPM_a)','STPM_b','std(STPM_b)'};
hdrs13 = {'t','std(t)','N','sampling (s)'};
hdrs14 = {'dt','mean count','std count','mean fit','std fit','GT'};
hdrs22 = {'t','std(t)'};
hdrs21 = {'GT states','GT TPM','states','TPM_w','-err_w','+err_w',...
    'TPM_wo','-err_wo','+err_wo','std(TPM_w)','std(-err_w)',...
    'std(+err_w)','std(TPM_wo)','std(-err_wo)','std(+err_wo)'};
hdrs1 = {'BIC','STPM','t_comp','dthista','dthistb'};
hdrs2 = {'TPM','t_comp'};
hdrs = {'DPH','BW'};

% show process
if sumexp
    disp('>> collect data for sum exp. fit model...')
else
    disp('>> collect data for DPH fit model...');
end

% read data files
dat = readData([rootdir,datdir],'model',sumexp);

% structure data for .json format
dat11 = [hdrs11; mat2rowcell(dat{2})];
dat12 = [hdrs12; dat{3}];
dat13 = [hdrs13; dat{4}];
dat14 = [hdrs14; mat2rowcell(dat{5}{1})];
dat15 = [hdrs14; mat2rowcell(dat{5}{2})];
dat1 = [hdrs1;{dat11,dat12,dat13,dat14,dat15}];
dat21 = [hdrs21;dat{1}];
dat22 = [hdrs22;dat{6}];
dat2 = [hdrs2;{dat21,dat22}];
dat = {dat1,dat2};

% write file
fname = 'data_figure3B';
if sumexp
    fname = [fname,'_sumexp'];
end
exportJson([destdir,'figure3',filesep,fname,'.json'],[hdrs;dat]);


function dat = readData(datdir,varargin)
% dat = readD(datdir)
% dat = readD(datdir,type)
%
% Gather analysis results through files and return them in an array.
%
% datdir: source directory
% type: (opt) 'wab' to read GT Nab/Nbb, 
%             'tau' to read GT state lifetimes
% dat: [nDat-by-5 to 9] data array

% defaults
ext_data = '_data.mat';
ext_simprm = '_simprm.mat';
ext_simres = '_simres.mat';
ext_dph = '_mldphres.mat';
ext_dphplot1 = '_dph_state1D*_dphplot';
ext_dphplot2 = '_dph_state2D*_dphplot';
ext_bw = '_bwres.mat';
type = 'none';
dtbinsim = 1;
dtbinexp = 10;

% collect read type
sumexp = false;
if ~isempty(varargin)
    type = varargin{1};
    if numel(varargin)>=2
        sumexp = varargin{2};
    end
end

if strcmp(type,'traj')
    fledat = importdata(datdir,'\t',3);
    dat = fledat.data(:,[1,5:8]);
    return
end

% check existence of simulated data
simlist = dir([datdir,filesep,'*',ext_simprm]);
issim = true;
if isempty(simlist)
	[~,subdir,~] = fileparts(datdir);
    simlist = dir([datdir,filesep,'*',ext_data]);
    if isempty(simlist)
        disp(['Neither simulation parameters or experimental data were ',...
            'found in sub-directory',subdir])
    else
        issim = false;
    end
    if ~issim
        disp(['No simulation parameters found for sub-directory',subdir])
    end
end

% initialize output
switch type
    case 'model'
        dat = cell(1,5);
        dat{5} = cell(1,2);
    case 'schm'
        dat = cell(1,3);
        dat{3} = {};
    otherwise
        dat = cell(1,2);
end

% read results
for f = 1:size(simlist)
    if issim
        % collect GT
        simprm = load([datdir,filesep,simlist(f,1).name]);
        val_b0 = simprm.val(end);
        val_a0 = simprm.val(1);

        % determine GT state degeneracies
        Db0 = sum(simprm.val==val_b0);
        Da0 = sum(simprm.val==val_a0);
        J0 = Da0+Db0;

        % calculates b-states exit probabilities
        wba0 = simprm.tp((Da0+1):end,1)./sum(simprm.tp((Da0+1):end,:),2);

        % calculates state lifetimes
        tau0 = (1./sum(simprm.tp,2))';
        taub0 = tau0(simprm.val==val_b0);

        % sort states according to lifetimes
        id0 = sortStates(simprm.val,tau0);
        ip0 = simprm.ip(id0);
        k0 = reorderMat(simprm.tp,id0);

        % calculates GT DPH parameters
        tp0 = k0;
        tp0(~~eye(J0)) = 1-sum(k0,2);
        [a0,T0,~] = calcDPHprm(ip0,tp0,Da0,Db0);

        % determines GT transition scheme for b states
        isip0_b = ~~(ip0((Da0+1):end)+sum(k0(1:Da0,(Da0+1):end),1));
        istp0_b = ~~k0((Da0+1):end,(Da0+1):end);
        isep0_b = ~~sum(k0((Da0+1):end,1:Da0),2);
        schmb0 = [false,isip0_b,false;[false(Db0,1),istp0_b,isep0_b];...
            false(1,Db0+2)];
        schm0 = k0;
        schm0(~~eye(size(schm0))) = 0;
        schm0 = double(schm0>0);
        
        % collects sample size and sampling rate
        N = simprm.N;
        rate = simprm.rate;
        dtbin = dtbinsim;
    else
        % collect experimental parameters
        expdat = load([datdir,filesep,simlist(f,1).name]);
        N = numel(expdat.dat.dt_obs);
        rate = expdat.prm.rate;
        dtbin = dtbinexp;
    end
	
	% import ML-DPH results
    if issim
        setname = simlist(f,1).name(1:end-length(ext_simprm));
        resfolder = [datdir,filesep,setname];
    else
        setname = simlist(f,1).name(1:end-length(ext_data));
        resfolder = datdir;
    end
    dphfle = dir([resfolder,filesep,setname,'_*',ext_dph]);
    if isempty(dphfle)
		disp(['No ML-DPH results found for set ',setname])
		continue
    end
    R = size(dphfle,1);
    
    % collect relevant files
    plotfle1 = dir([resfolder,filesep,setname,'_*',ext_dphplot1]);
    plotfle2 = dir([resfolder,filesep,setname,'_*',ext_dphplot2]);
    if issim
        datfle = dir([resfolder,filesep,setname,'_*',ext_simres]);
    else
        datfle = dir([resfolder,filesep,setname,ext_data]);
    end

	% initializes tables
	Db = [];
	acc_dph = [];
	acc_bw = [];
    TPMs = [];
    BICs = [];
    t_dph = [];
    t_bw = [];
    STPMs = cell(1,2);
    if strcmp(type,'model')
        plotdat = cell(1,2);
        success = 0;
    end
    
    % collect results for each replicate
	for r = 1:R
		dphres = load([dphfle(r,1).folder,filesep,dphfle(r,1).name]);
        if sumexp
            dphres = dphres.expres;
        else
            dphres = dphres.dphres;
        end
        tp_a = dphres{3}.tp_fit{1};
        k_a = tp_a;
        k_a(~~eye(size(k_a))) = 0;
        tp_b = dphres{3}.tp_fit{2};
        k_b = tp_b;
        k_b(~~eye(size(k_b))) = 0;
        schmb = dphres{3}.schm{2};
        
        val = sort(unique(dphres{4}));
        val_a = val(1);
        val_b = val(2);
        if isempty(k_b) % ML-DPH did not converge
            continue
        end
        
		Db = cat(2,Db,size(tp_b,1));
		Da = size(tp_a,1);
		
        if ~strcmp(type,'model')
            % calculates model fidelity after ML-DPH
            if Db(end)~=Db0
                acc_r = 0;
            else
                acc_r = [];
                for db = 1:Db
                    idb = circshift(1:Db,db);
                    isipb = schmb(1,2:end-1);
                    isepb = schmb(2:end-1,end);
                    schmschft = schmb;
                    schmschft(1,2:end-1) = isipb(idb);
                    schmschft(2:end-1,end) = isepb(idb);
                    schmschft(2:end-1,2:end-1) = ...
                        reorderMat(schmb(2:end-1,2:end-1),idb);
                    acc_r = cat(2,acc_r,calcmodelacc(schmb0,schmschft));
                end
            end
            acc_dph = cat(2,acc_dph,max(acc_r));
            
            % calculates model fidelity after BW
            if strcmp(type,'schm')
                bwfle = [resfolder,filesep,setname,'_',num2str(r),ext_bw];
                if ~exist(bwfle,'file')
                    disp(['No BW results found for set ',setname])
                    continue
                end
                acc_r = [];
                if ~(Db(end)==Db0 && Da==Da0)
                    acc_r = 0;
                else
                    bwres = load(bwfle);
                    tp_wo = bwres.bwres_wo{1};
                    tp_wo(~~eye(size(tp_wo))) = 0;
                    tauwo = 1./sum(tp_wo,2);
                    states = dphres{4};
                    id_wo = sortStates(states,tauwo);
                    tp_wo = reorderMat(tp_wo,id_wo);
                    schm_wo = double(tp_wo>0);
                    for db = 1:Db
                        id = [1:Da,circshift(Da+1:Da+Db,db)];
                        schmschft = reorderMat(schm_wo,id);
                        acc_r = cat(2,acc_r,calcmodelacc(schm0,schmschft));
                    end
                end
                acc_bw = cat(2,acc_bw,max(max(acc_r)));
            end
            
            if r==1
                simres = ...
                    load([resfolder,filesep,datfle(r,1).name]);
                dt_gt = simres.res.dt_gt;
                for n = 1:numel(dt_gt)
                    dt_gt{n}(:,1) = dt_gt{n}(:,1)*simprm.rate;
                end
                plotdat = bindthist(importdata([resfolder,filesep,...
                    plotfle2(r,1).name]),dtbin);
            end
            
        elseif ~issim || (Db(end)==Db0 && Da==Da0) % only for success
            % get computation time
            t_dph = cat(2,t_dph,dphres{3}.t_dphtest);
            
            % collect ML-DPH sub-transition matrices (STPMs)
            tau_a = (1./sum(k_a,2))';
            tau_b = (1./sum(k_b,2))';
            STPMs{1} = cat(3,STPMs{1},reorderMat(tp_a,...
                sortStates(repmat(val_a,1,Da),tau_a)));
            STPMs{2} = cat(3,STPMs{2},reorderMat(tp_b,...
                sortStates(repmat(val_b,1,Db(end)),tau_b)));
            
            % collect BW transition probability matrices (TPMs)
            bwfle = [resfolder,filesep,setname,'_',num2str(r),ext_bw];
            if ~exist(bwfle,'file')
                disp(['No BW results found for set ',setname])
                continue
            end
            bwres = load(bwfle);
            states = dphres{4};
            
            % get computation time
            if numel(bwres.bwres_wo)>=5
                t_bw = cat(2,t_bw,bwres.bwres_wo{5});
            end
            
            % use same state order (ascending value, then ascending lifetime)
            tp_wo = bwres.bwres_wo{1};
            tp_wo(~~eye(size(tp_wo))) = 0;
            tauwo = 1./sum(tp_wo,2);
            id_wo = sortStates(states,tauwo);
            if issim
                tp_w = bwres.bwres_w{1};
                tp_w(~~eye(size(tp_w))) = 0;
                tauw = 1./sum(tp_w,2);
                id_w = sortStates(states,tauw);
                TPMs = cat(3,TPMs,[simprm.val(id0)',tp0,states(id_w)',...
                    reorderMat(bwres.bwres_w{1},id_w),...
                    reorderMat(bwres.bwres_w{2}(:,:,1),id_w),...
                    reorderMat(bwres.bwres_w{2}(:,:,2),id_w),...
                    reorderMat(bwres.bwres_wo{1},id_wo),...
                    reorderMat(bwres.bwres_wo{2}(:,:,1),id_wo), ...
                    reorderMat(bwres.bwres_wo{2}(:,:,2),id_wo)]);
            else
                J = numel(states);
                TPMs = cat(3,TPMs,[states(id_wo)',...
                    reorderMat(bwres.bwres_wo{1},id_wo),...
                    reorderMat(bwres.bwres_wo{2}(:,:,1),id_wo), ...
                    reorderMat(bwres.bwres_wo{2}(:,:,2),id_wo)]);
            end
            fitmdl = dphres{2}(dphres{2}(:,1)==1,[2,3]);
            bicdat = [];
            for s = 1:size(fitmdl,1)
                ismdl_s = dphres{2}(:,2)==fitmdl(s,1) & ...
                    dphres{2}(:,3)==fitmdl(s,2);
                if any(dphres{2}(dphres{2}(:,1)==1 & ismdl_s,6)==0)
                    bic_a = Inf;
                else
                    bic_a = dphres{2}(dphres{2}(:,1)==1 & ismdl_s,5);
                end
                if any(dphres{2}(dphres{2}(:,1)==2 & ismdl_s,6)==0)
                    bic_b = Inf;
                else
                    bic_b = dphres{2}(dphres{2}(:,1)==2 & ismdl_s,5);
                end
                bicdat = cat(1,bicdat,[fitmdl(s,[1,2]),bic_a,bic_b]);
            end
            BICs = cat(3,BICs,bicdat);

            dthista = importdata([resfolder,filesep,plotfle1(r,1).name]);
            dthistb = importdata([resfolder,filesep,plotfle2(r,1).name]);
            plotdat{1} = cat(1,plotdat{1},...
                [dthista,r*ones(size(dthista,1),1)]);
            plotdat{2} = cat(1,plotdat{2},...
                [dthistb,r*ones(size(dthistb,1),1)]);
            success = success+1;
        end
	end
	
	% calculate std
    switch type
        case 'tau'
            dat{1} = cat(1,dat{1},[Db0,sort(taub0),mean(Db),std(Db),...
                mean(acc_dph(acc_dph>0)),std(acc_dph(acc_dph>0))]);
            dthist = dthiststats(sortrows(plotdat,1),dtbin,1);
            dat{2} = cat(1,dat{2},{[dthist,...
                calcDPHprob(T0{2},a0{2},dthist(:,1))]});
            
        case 'tauwba'
            fact = sort(taub0);
            fact = round(fact(2)/fact(1));
            dat{1} = cat(1,dat{1},[Db0,fact,sort(taub0),mean(wba0),...
                mean(Db),std(Db)]);
            dthist = dthiststats(sortrows(plotdat,1),dtbin,1);
            dat{2} = cat(1,dat{2},{[dthist,...
                calcDPHprob(T0{2},a0{2},dthist(:,1))]});
            
        case 'wba'
            dat{1} = cat(1,dat{1},[Db0,mean(wba0),mean(Db),std(Db),...
                mean(acc_dph(acc_dph>0)),std(acc_dph(acc_dph>0))]);
            dthist = dthiststats(sortrows(plotdat,1),dtbin,1);
            dat{2} = cat(1,dat{2},{[dthist,...
                calcDPHprob(T0{2},a0{2},dthist(:,1))]});
        case 'schm'
            dat{1} = cat(1,dat{1},[Db0,mean(Db),std(Db),mean(acc_bw(acc_bw>0)),...
                std(acc_bw(acc_bw>0))]);
            dat{2} = cat(1,dat{2},schm0);
            dthist = dthiststats(sortrows(plotdat,1),dtbin,1);
            dat{3} = cat(1,dat{3},{[dthist,...
                calcDPHprob(T0{2},a0{2},dthist(:,1))]});
        case 'model'
            TPMs_std = std(TPMs,[],3);
            TPMs = mean(TPMs,3);
            STPMs = {mean(STPMs{1},3),std(STPMs{1},[],3),...
                mean(STPMs{2},3),std(STPMs{2},[],3)};
            bica = permute(BICs(:,3,:),[1,3,2]);
            bicb = permute(BICs(:,4,:),[1,3,2]);
            S = size(BICs,1);
            bicstats = zeros(S,4);
            for s = 1:S
                bicstats(s,1) = mean(bica(s,~isinf(bica(s,:))));
                bicstats(s,2) = mean(bicb(s,~isinf(bicb(s,:))));
                bicstats(s,3) = std(bica(s,~isinf(bica(s,:))));
                bicstats(s,4) = std(bicb(s,~isinf(bicb(s,:))));
            end
            BICs = [BICs(:,[1,2],1),bicstats];
            if issim
                dat{1} = {TPMs(:,1),TPMs(:,2:J0+1),TPMs(:,J0+2),...
                    TPMs(:,(J0+3):(2*J0+2)),TPMs(:,(2*J0+3):(3*J0+2)),...
                    TPMs(:,(3*J0+3):(4*J0+2)),TPMs(:,(4*J0+3):(5*J0+2)),...
                    TPMs(:,(5*J0+3):(6*J0+2)),TPMs(:,(6*J0+3):(7*J0+2)),...
                    TPMs_std(:,(J0+3):(2*J0+2)),TPMs_std(:,(2*J0+3):(3*J0+2)),...
                    TPMs_std(:,(3*J0+3):(4*J0+2)),TPMs_std(:,(4*J0+3):(5*J0+2)),...
                    TPMs_std(:,(5*J0+3):(6*J0+2)),TPMs_std(:,(6*J0+3):(7*J0+2))};
            else
                dat{1} = {TPMs(:,1),TPMs(:,2:(J+1)),TPMs(:,(J+2):(2*J+1)),...
                    TPMs(:,(2*J+2):(3*J+1)),TPMs_std(:,2:(J+1)),...
                    TPMs_std(:,(J+2):(2*J+1)),TPMs_std(:,(2*J+2):(3*J+1))};
            end
            dat{2} = BICs;
            dat{3} = STPMs;
            dat{4} = {mean(t_dph),std(t_dph),N,1/rate};
            dthista = dthiststats(bindthist(plotdat{1},dtbin),dtbin,...
                success);
            dthistb = dthiststats(bindthist(plotdat{2},dtbin),dtbin,...
                success);
            if issim
                dat{5}{1} = [dthista,calcDPHprob(T0{1},a0{1},dthista(:,1))];
                dat{5}{2} = [dthistb,calcDPHprob(T0{2},a0{2},dthistb(:,1))];
            else
                dat{5}{1} = dthista;
                dat{5}{2} = dthistb;
            end
            dat{6} = {mean(t_bw),std(t_bw)};
            return
            
        case 'none'
            dat{1} = cat(1,dat{1},...
                [Db0,mean(Db),std(Db),mean(acc_dph(acc_dph>0)),std(acc_dph(acc_dph>0))]);
            dthist = dthiststats(sortrows(plotdat,1),dtbin,1);
            dat{2} = cat(1,dat{2},{[dthist,...
                calcDPHprob(T0{2},a0{2},dthist(:,1))]});
    end
end
switch type
    case 'none'
        if ~isempty(dat{1})
            [~,id] = sort(dat{1}(:,1));
            dat{1} = dat{1}(id,:);
            dat{2} = dat{2}(id,:);
        end
    case 'tau'
        if ~isempty(dat{1})
            [~,id] = sort(dat{1}(:,3));
            dat{1} = dat{1}(id,:);
            dat{2} = dat{2}(id,:);
        end
    case 'wba'
        if ~isempty(dat{1})
            [~,id] = sort(dat{1}(:,2));
            dat{1} = dat{1}(id,:);
            dat{2} = dat{2}(id,:);
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

[R,C] = size(mat);
rowcell = mat2cell(mat,R,ones(1,C));



function acc = calcmodelacc(schm0,schm)
% acc = calcmodelacc(schm0,schm)
%
% Calculates accuracy of transition scheme
%
% schm0: [D+2-by-D+2] GT transition scheme
% schm: [D+2-by-D+2] inferred transition scheme
% acc: model accuracy

TP = sum(sum(schm0 & schm));
TN = sum(sum(~schm0 & ~schm));
FP = sum(sum(~schm0 & schm));
FN = sum(sum(schm0 & ~schm));
acc = (TP+TN)/(TP+TN+FP+FN);


function st = dthiststats(dthist,dtbin,R)
% st = dthiststats(dthist,dtbin,R)
%
% Calculates mean and standard deviaition of dwell time counts and DPH fit
% among simulation replicates
%
% dthist: [R*nbin-by-3] dwell time, counts, fit
% dtbin: bin size
% R: number of replicates
% st: [nbin-by-5] dwell time, mean counts, std counts, mean fit, std fit

% defaults
ngtdat = 50; % number of GT PDF data points

st = [];
maxdt = ceil(max(dthist(:,1))/dtbin);

% bins to calculate GT PDF for even if dwell time count is null
gtbin = unique(round(logspace(log10(1),log10(maxdt),ngtdat)));

for bin = 1:maxdt
    id = find(dthist(:,1)==bin*dtbin);
    if isempty(id) && bin>1
        continue
    end
    cnt = dthist(id,2)';
    if sum(cnt)==0 && all(bin~=gtbin)
        continue
    end
    prob = dthist(id,3)';
    if numel(id)<R
        cnt = cat(2,cnt,zeros(1,R-numel(id)));
        prob = cat(2,prob,zeros(1,R-numel(id)));
    end
    st = cat(1,st,[bin*dtbin,mean(cnt),std(cnt),mean(prob),...
        std(prob)]);
end


function [a0,T0,t0] = calcDPHprm(ip,tp,Da,Db)
% [a0,T0,t0] = calcDPHprm(ip,tp,Da,Db)
% 
% Calculates theoretical parameters of discrete phase-type distributions
% of two observed states.
%
% ip: [1-by-Da+Db] initial state probabilties
% tp: [Da+Db-by-Da+Db] transition probabilities
% Da: state degeneracy of observed state 1
% Db: state degeneracy of observed state 2
% a0: {1-by-2}[1-by-D] initial state probabilities
% T0: {1-by-2}[D-by-D] generator
% t0: {1-by-2}[D-by-1] exit probabilties

% initializes output
a0 = cell(1,2);
T0 = cell(1,2);

T0{1} = tp(1:Da,1:Da);
t0{1} = sum(tp(1:Da,(Da+1):end),2);
a0{1} = mean(tp((Da+1):end,1:Da)./...
    repmat(sum(tp((Da+1):end,1:Da),2),[1,Da]),1);
a0{1} = a0{1}/sum(a0{1});

T0{2} = tp((Da+1):end,(Da+1):end);
t0{2} = sum(tp((Da+1):end,1:Da),2);
a0{2} = mean(tp(1:Da,(Da+1):end)./...
    repmat(sum(tp(1:Da,(Da+1):end),2),[1,Db]),1);
a0{2} = a0{2}/sum(a0{2});


function [a0,T0] = calcexpDPHprm(dt,Da,Db)
% [a0,T0] = calcexpDPHprm(dt,Da,Db)
%
% Calculates experimental parameters of discrete phase-type distributions
% of two observed states.
%
% dt: [ndt-by-3] dwell time, state index, index of state after transition
% Da: state degeneracy of observed state 1
% Db: state degeneracy of observed state 2
% a0: {1-by-2}[1-by-D] experimental initial state probabilities
% T0: {1-by-2}[D-by-D] experimental generators

if ~iscell(dt)
    dt = {dt};
end
a0{1} = zeros(1,Da);
a0{2} = zeros(1,Db);
T0{1} = zeros(Da,Da+1);
T0{2} = zeros(Db,Db+1);
tcum{1} = zeros(Da,1);
tcum{2} = zeros(Db,1);
for n = 1:numel(dt)
    id_a = find(dt{n}(:,2)>=1 & dt{n}(:,2)<=Da);
    if ~isempty(id_a)
        frsta = dt{n}(id_a(1),2);
        if numel(id_a)>1
            frsta = cat(1,frsta,dt{n}(id_a(find(diff(id_a)>1)+1),2));
        end
    else
        frsta = [];
    end
    id_b = find(dt{n}(:,2)>=(Da+1) & dt{n}(:,2)<=(Da+Db));
    if ~isempty(id_b)
        frstb = dt{n}(id_b(1),2);
        if numel(id_b)>1
            frstb = cat(1,frstb,dt{n}(id_b(find(diff(id_b)>1)+1),2));
        end
    else
        frstb = [];
    end
    for a = 1:Da
        a0{1}(a) = a0{1}(a)+sum(frsta==a);
    end
    for b = 1:Db
        a0{2}(b) = a0{2}(b)+sum(frstb==(Da+b));
    end
    for a1 = 1:Da
        for a2 = 1:Da
            if a1==a2
                continue
            end
            T0{1}(a1,a2) = T0{1}(a1,a2) + ...
                sum(dt{n}(:,2)==(a1) & dt{n}(:,3)==(a2));
        end
        T0{1}(a1,end) = T0{1}(a1,end) + ...
            sum(dt{n}(:,2)==(a1) & dt{n}(:,3)>Da);
        tcum{1}(a1,1) = tcum{1}(a1,1)+sum(dt{n}(dt{n}(:,2)==a1,1));
    end
    for b1 = 1:Db
        for b2 = 1:Db
            if b1==b2
                continue
            end
            T0{2}(b1,b2) = T0{2}(b1,b2) + ...
                sum(dt{n}(:,2)==(Da+b1) & dt{n}(:,3)==(Da+b2));
        end
        T0{2}(b1,end) = T0{2}(b1,end) + ...
            sum(dt{n}(:,2)==(Da+b1) & dt{n}(:,3)<=Da);
        tcum{2}(b1,1) = tcum{2}(b1,1)+sum(dt{n}(dt{n}(:,2)==(Da+b1),1));
    end
end
a0{1} = a0{1}/sum(a0{1});
a0{2} = a0{2}/sum(a0{2});
T0{1} = T0{1}./repmat(tcum{1},1,Da+1);
T0{2} = T0{2}./repmat(tcum{2},1,Db+1);
T0{1}(~~eye(size(T0{1}))) = 1-sum(T0{1},2);
T0{2}(~~eye(size(T0{2}))) = 1-sum(T0{2},2);
T0{1} = T0{1}(:,1:end-1);
T0{2} = T0{2}(:,1:end-1);


function prob = calcDPHprob(T,ip,dt)
% prob = calcDPHprob(T,ip,dt)
%
% Calculates probability density of discrete phase-type distribution from
% input parameters.
% 
% T: [D-by-D] generator
% ip: [1-by-D] initial state probability
% dt: vector of times to calculate probability densities for
% prob: vector of probability densities

prob = zeros(size(dt));
D = size(T,1);
v_e = ones(D,1);
t = v_e-T*v_e;
for n = 1:numel(dt)
    prob(n) = ip*(T^(dt(n)-1))*t;
end


function dthist = bindthist(dthist0,dtbin)
% dthist = bindthist(dthist0,dtbin,T,ip)
% 
% Over-bin dwell time histogram counts.
%
% dthist0: [nbin0-by-3 or 4] dwell time histogram (times, counts, DPH fit,
%   replicate)
% dtbin: how many original bins must be over-binned
% dthist: [nbin-by-3 or 4] over-binned dwell time histogram (times, counts, 
%   DPH fit, replicate)

% initializes output
dthist = [];

% get replicate data location
ndt = size(dthist0,1);
if size(dthist0,2)>=4
    rs = unique(dthist0(:,4));
    R = numel(rs);
    id_r = false(ndt,R);
    for r = 1:R
        id_r(:,r) = dthist0(:,4)==rs(r);
    end
else
    R = 1;
    id_r = true(ndt,1);
end

% over-bins histogram
maxdt = max(dthist0(:,1));
for bin = 1:ceil(maxdt/dtbin)
    if bin==1
        lowlim = 0;
    else
        lowlim = (bin-1)*dtbin;
    end
    id = find(dthist0(:,1)>lowlim & dthist0(:,1)<=bin*dtbin);
    if isempty(id)
        continue
    end
    for r = 1:R
        dthist = cat(1,dthist,[bin*dtbin,sum(dthist0(id(id_r(id,r)),2)),...
            sum(dthist0(id(id_r(id,r)),3))]);
    end
end


function mat = reorderMat(mat0,id)
% mat = reorderMat(mat0,id)
%
% Order rows and columns of a square matrix according to input order.
%
% mat0: [N-by-N] original matrix
% id: [1-by-N] rows and columns new order
% mat: [N-by-N] ordered matrix

J = size(mat0,1);
if numel(id)~=J
    disp(['reorderMat: the number of indexes is different than the matrix',...
        ' size.'])
    mat = [];
    return
end
mat = zeros(size(mat0));
for j1 = 1:J
    for j2 = 1:size(mat0,2)
        if j2>J
            mat(j1,j2) = mat0(id(j1),j2);
        else
            mat(j1,j2) = mat0(id(j1),id(j2));
        end
    end
end


function id = sortStates(stateval,tau)
% id = sortStates(stateval,tau)
%
% Sort input degenerate states according to their lifetimes.
% Sorting is first performed on state values, then on lifetimes.
%
% stateval: [1-by-J] state values
% tau: [1-by-J] state lifetimes
% id: [1-by-J] new state order

id = [];
val = sort(unique(stateval));
for v = 1:numel(val)
    idv = find(stateval==val(v));
    [~,idtau] = sort(tau(idv));
    id = cat(2,id,idv(idtau));
end


function exportJson(jfile,dat)
% exportJson(jfile,dat)
%
% Export data into a .json file using JSON formatting
%
% jfile: destination file address
% dat: {2-by-C} header strings and numerical data columns

[pname,~,~] = fileparts(jfile);
if ~exist(pname,'dir')
    mkdir(pname);
end

f = fopen(jfile,'w');
writeJson(f,dat,0);
fclose(f);
disp(['File ',jfile,' was successfully exported!']);


function writeJson(f,dat,ntab)
% writeJson(f,dat,ntab)
%
% Write header and data into an ASCII file using JSON formatting
%
% f: file identifier
% dat: {2-by-C} header strings (1st row) and associated numerical data (2nd 
%  row) or [R-by-C] numerical data.
% ntab: data indentation

C = size(dat,2);
tabs = repmat('\t',[1,ntab]);
for c = 1:C
    if c==1 
        if size(dat,1)==2
            fprintf(f,'{');
        else
            fprintf(f,'[');
        end
    end
    fprintf(f,['\n',tabs,'\t']);
    if size(dat,1)==2
        if ~isempty(dat{1,c}) && ischar(dat{1,c}) % write header
            fprintf(f,['\"',dat{1,c},'\": ']);
        end
        if iscell(dat{2,c})
            writeJson(f,dat{2,c},ntab+1);
        else
            [R,ndat] = size(dat{2,c});
            if R>1
                fprintf(f,'[');
            end
            for r = 1:R
                if ndat>0
                    fmtdat = '%d';
                    if ndat>1
                        fprintf(f,'[');
                        fmtdat = cat(2,fmtdat,repmat(',%d',[1,ndat-1]));
                    end
                    fprintf(f,fmtdat,dat{2,c}(r,:));
                    if ndat>1
                        fprintf(f,']');
                    end
                end
                if r~=R
                    fprintf(f,',');
                end
            end
            if R>1
                fprintf(f,']');
            end
        end
    else
        [R,ndat] = size(dat);
        if iscell(dat)
            for r = 1:R
                writeJson(f,dat{r,1},ntab+1);
                if r~=R
                    fprintf(f,',');
                end
            end
        else
            if R>1
                fprintf(f,'[');
            end
            for r = 1:R
                if ndat>0
                    fmtdat = '%d';
                    if ndat>1
                        fprintf(f,'[');
                        fmtdat = cat(2,fmtdat,repmat(',%d',[1,ndat-1]));
                    end
                    fprintf(f,fmtdat,dat(r,:));
                    if ndat>1
                        fprintf(f,']');
                    end
                end
                if r~=R
                    fprintf(f,',');
                end
            end
            if R>1
                fprintf(f,']');
            end
        end
    end
    if c~=C
        fprintf(f,',');
    end
    if c==C
        if size(dat,1)==2
            fprintf(f,['\n',tabs,'}']);
        else
            fprintf(f,['\n',tabs,']']);
        end
    end
end

