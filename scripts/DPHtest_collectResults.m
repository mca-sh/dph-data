function DPHtest_collectResults(rootdir,destdir,varargin)
% DPHtest_collectResults(rootdir,destdir)
% DPHtest_collectResults(rootdir,destdir,subdir)
%
% Collect results for plotting figures and epxort them to .json file format.
%
% rootdir: source directory
% destdir: destination directory
% subdir: cell array containing specific subdirectories to export

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
else
    subfig = {'figure2AB','figure2C','figure2D','figure3A','figure3B',...
        'figureS1A'};
end
strsubfig = '';
for sd = 1:numel(subfig)
    strsubfig = cat(2,strsubfig,subfig{sd},' ');
end
disp(['files will be exported for: ',strsubfig(1:end-1),'.']);

% figure 2AB: D and model fidelity = f(D_GT) for models I and II
if any(contains(subfig,'figure2AB'))
    disp('process subfigure 2AB...')
    datdirI = 'dataset1';
    datdirII = 'dataset2';
    hdrs_2AB = {'GT Db_I','Db_I','std(Db_I)','MF_I','std(MF_I)','GT Db_II',...
        'Db_II','std(Db_II)','MF_II','std(MF_II)'};
    
    disp('>> collect data for DPH fit model...')
    dat_2AB = [readData([rootdir,datdirI]),readData([rootdir,datdirII])];
    dat_2AB = mat2cell(dat_2AB,size(dat_2AB,1),ones(1,size(dat_2AB,2)));
    exportJson([destdir,'figure2',filesep,'data_figure2AB.json'],...
        [hdrs_2AB;dat_2AB]);
    
    disp('>> collect data for sum exp. fit model...')
    dat_2AB = [readData([rootdir,datdirI],'none',true),...
        readData([rootdir,datdirII],'none',true)];
    dat_2AB = mat2cell(dat_2AB,size(dat_2AB,1),ones(1,size(dat_2AB,2)));
    exportJson([destdir,'figure2',filesep,'data_figure2AB_sumexp.json'],...
        [hdrs_2AB;dat_2AB]);
end

% figure 2C: D and model fidelity = f(Nab/Nbb) for model II and three DPH fit curves
if any(contains(subfig,'figure2C'))
    disp('process subfigure 2C...')
    datdir = 'dataset3';
    hdrs_2C2 = {'dt','summed count','std count','summed fit','std fit',...
        'summed mdl','std mdl'};
    hdrs_2C = {'GT Db_II','Nab/Nbb','Db_II','std(Db_II)','MF_II',...
        'std(MF_II)','dthist'};
    
    disp('>> collect data for DPH fit model...')
    dat = readData([rootdir,datdir],'wba');
    dat_2C1 = mat2cell(dat{1},size(dat{1},1),ones(1,size(dat{1},2)));
    dat_2C2 = {};
    for n = 1:numel(dat{2})
        dat_2C2 = cat(1,dat_2C2,{[hdrs_2C2;mat2cell(dat{2}{n},...
            size(dat{2}{n},1),ones(1,size(dat{2}{n},2)))]});
    end
    dat_2C = [dat_2C1,{dat_2C2}];
    exportJson([destdir,'figure2',filesep,'data_figure2C.json'],...
        [hdrs_2C;dat_2C]);
    
    disp('>> collect data for sum exp. fit model...')
    dat = readData([rootdir,datdir],'wba',true);
    dat_2C1 = mat2cell(dat{1},size(dat{1},1),ones(1,size(dat{1},2)));
    dat_2C2 = {};
    for n = 1:numel(dat{2})
        dat_2C2 = cat(1,dat_2C2,{[hdrs_2C2;mat2cell(dat{2}{n},...
            size(dat{2}{n},1),ones(1,size(dat{2}{n},2)))]});
    end
    dat_2C = [dat_2C1,{dat_2C2}];
    exportJson([destdir,'figure2',filesep,'data_figure2C_sumexp.json'],...
        [hdrs_2C;dat_2C]);
end

% figure 2D: D and model fidelity = f(tau_b2) for model I and II
if any(contains(subfig,'figure2D'))
    disp('process subfigure 2D...')
    datdirI = 'dataset4';
    datdirII = 'dataset5';
    hdrs_2D = {'GT Db_I','GT tau_b1_I','GT tau_b2_I','GT tau_b3_I','Db_I',...
        'std(Db_I)','MF_I','std(MF_I)','GT Db_II','GT tau_b1_II',...
        'GT tau_b2_II','GT tau_b3_II','Db_II','std(Db_II)','MF_II',...
        'std(MF_II)'};
    
    disp('>> collect data for DPH fit model...')
    dat_2D = [readData([rootdir,datdirI],'tau'),...
        readData([rootdir,datdirII],'tau')];
    dat_2D = mat2cell(dat_2D,size(dat_2D,1),ones(1,size(dat_2D,2)));
    exportJson([destdir,'figure2',filesep,'data_figure2D.json'],...
        [hdrs_2D;dat_2D]);
    
    disp('>> collect data for sum exp. fit model...')
    dat_2D = [readData([rootdir,datdirI],'tau',true),...
        readData([rootdir,datdirII],'tau',true)];
    dat_2D = mat2cell(dat_2D,size(dat_2D,1),ones(1,size(dat_2D,2)));
    exportJson([destdir,'figure2',filesep,'data_figure2D_sumexp.json'],...
        [hdrs_2D;dat_2D]);
end

% figure 3A: GT and inferred TPMs for quenched disorder, ML-DPH BICs and TPs, and computation time
if any(contains(subfig,'figure3A'))
    disp('process subfigure 3A...')
    datdir = 'dataset6';
    hdrs_3A11 = {'D','scheme','BIC_a','BIC_b','std(BIC_a)','std(BIC_b)'};
    hdrs_3A12 = {'STPM_a','std(STPM_a)','STPM_b','std(STPM_b)'};
    hdrs_3A13 = {'t','std(t)'};
    hdrs_3A1 = {'BIC','STPM','t_comp'};
    hdrs_3A2 = {'GT states','GT TPM','states','TPM_w','-err_w','+err_w',...
        'TPM_wo','-err_wo','+err_wo','std(TPM_w)','std(-err_w)',...
        'std(+err_w)','std(TPM_wo)','std(-err_wo)','std(+err_wo)'};
    hdrs_3A = {'DPH','BW'};
    
    disp('>> collect data for DPH fit model...')
    dat = readData([rootdir,datdir],'model');
    dat_3A11 = [hdrs_3A11; ...
        mat2cell(dat{2},size(dat{2},1),ones(1,size(dat{2},2)))];
    dat_3A12 = [hdrs_3A12; dat{3}];
    dat_3A13 = [hdrs_3A13; dat{4}];
    dat_3A1 = [hdrs_3A1;{dat_3A11,dat_3A12,dat_3A13}];
    dat_3A2 = [hdrs_3A2;dat{1}];
    dat_3A = {dat_3A1,dat_3A2};
    exportJson([destdir,'figure3',filesep,'data_figure3A.json'],...
        [hdrs_3A;dat_3A]);
    
    disp('>> collect data for sum exp. fit model...')
    dat = readData([rootdir,datdir],'model',true);
    dat_3A11 = [hdrs_3A11; ...
        mat2cell(dat{2},size(dat{2},1),ones(1,size(dat{2},2)))];
    dat_3A12 = [hdrs_3A12; dat{3}];
    dat_3A13 = [hdrs_3A13; dat{4}];
    dat_3A1 = [hdrs_3A1;{dat_3A11,dat_3A12,dat_3A13}];
    dat_3A2 = [hdrs_3A2;dat{1}];
    dat_3A = {dat_3A1,dat_3A2};
    exportJson([destdir,'figure3',filesep,'data_figure3A_sumexp.json'],...
        [hdrs_3A;dat_3A]);
end


% figure 3B: GT and inferred TPM for dynamic disorder, ML-DPH BICs and TPs and computation time
if any(contains(subfig,'figure3B'))
    disp('process subfigure 3B...')
    datdir = 'dataset7';
    hdrs_3A11 = {'D','scheme','BIC_a','BIC_b','std(BIC_a)','std(BIC_b)'};
    hdrs_3A12 = {'STPM_a','std(STPM_a)','STPM_b','std(STPM_b)'};
    hdrs_3A13 = {'t','std(t)'};
    hdrs_3A1 = {'BIC','STPM','t_comp'};
    hdrs_3A2 = {'GT states','GT TPM','states','TPM_w','-err_w','+err_w',...
        'TPM_wo','-err_wo','+err_wo','std(TPM_w)','std(-err_w)',...
        'std(+err_w)','std(TPM_wo)','std(-err_wo)','std(+err_wo)'};
    hdrs_3A = {'DPH','BW'};
    
    disp('>> collect data for DPH fit model...')
    dat = readData([rootdir,datdir],'model');
    dat_3B11 = [hdrs_3A11; ...
        mat2cell(dat{2},size(dat{2},1),ones(1,size(dat{2},2)))];
    dat_3B12 = [hdrs_3A12; dat{3}];
    dat_3B13 = [hdrs_3A13; dat{4}];
    dat_3B1 = [hdrs_3A1;{dat_3B11,dat_3B12,dat_3B13}];
    dat_3B2 = [hdrs_3A2;dat{1}];
    dat_3B = {dat_3B1,dat_3B2};
    exportJson([destdir,'figure3',filesep,'data_figure3B.json'],...
        [hdrs_3A;dat_3B]);
    
    disp('>> collect data for sum exp. fit model...')
    dat = readData([rootdir,datdir],'model',true);
    dat_3B11 = [hdrs_3A11; ...
        mat2cell(dat{2},size(dat{2},1),ones(1,size(dat{2},2)))];
    dat_3B12 = [hdrs_3A12; dat{3}];
    dat_3B13 = [hdrs_3A13; dat{4}];
    dat_3B1 = [hdrs_3A1;{dat_3B11,dat_3B12,dat_3B13}];
    dat_3B2 = [hdrs_3A2;dat{1}];
    dat_3B = {dat_3B1,dat_3B2};
    exportJson([destdir,'figure3',filesep,'data_figure3B_sumexp.json'],...
        [hdrs_3A;dat_3B]);
end

% figure 4B: EBS-IBS intensity- and FRET-time traces
% figure 4C: EBS-IBS dwell time histograms, ML-DPH BICs and best fit curves
% figure 4D: EBS-IBS inferred TPM and computation time

disp('Process completed!');


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
ext_simprm = '_simprm.mat';
ext_simres = '_simres.mat';
ext_dph = '_mldphres.mat';
ext_dphplot = '_dph_state2D*_dphplot';
ext_bw = '_bwres.mat';
type = 'none';
dtbin = 10;

% collect read type
sumexp = false;
if ~isempty(varargin)
    type = varargin{1};
    if numel(varargin)>=2
        sumexp = varargin{2};
    end
end

% check existence of simulated data
simlist = dir([datdir,filesep,'*',ext_simprm]);
if isempty(simlist)
	[~,subdir,~] = fileparts(datdir);
	disp(['No simulation parameters found for sub-directory',subdir])
end

% initialize output
switch type
    case 'model'
        dat = cell(1,4);
    case 'wba'
        dat = cell(1,2);
    otherwise
        dat = [];
end

% read results
for f = 1:size(simlist)
	% collect GT
	simprm = load([datdir,filesep,simlist(f,1).name]);
	val_b = simprm.val(end);
	val_a = simprm.val(1);
	Db0 = sum(simprm.val==val_b);
    Da0 = sum(simprm.val==val_a);
	wba0 = simprm.tp(2:end,1)./sum(simprm.tp(2:end,:),2);
    tau0 = (1./sum(simprm.tp,2))';
    taub0 = tau0(simprm.val==val_b);
    id0 = sortStates(simprm.val,tau0);
    tp0 = reorderMat(simprm.tp,id0);
    J = numel(id0);
	
	% collect Db
	setname = simlist(f,1).name(1:end-length(ext_simprm));
	dphfle = ...
        dir([datdir,filesep,setname,filesep,setname,'_*',ext_dph]);
    if isempty(dphfle)
		disp(['No ML-DPH results found for set ',setname])
		continue
    end
    if strcmp(type,'wba')
        plotfle = dir([datdir,filesep,setname,filesep,setname,'_*',...
            ext_dphplot]);
        simfle = dir([datdir,filesep,setname,filesep,setname,'_*',...
            ext_simres]);
    end
	R = size(dphfle,1);
	Db = zeros(1,R);
	fdl = zeros(1,R);
    TPMs = [];
    BICs = [];
    t_comp = zeros(1,R);
    STPMs = cell(1,2);
    plotdat = [];
	for r = 1:R
		dphres = load([datdir,filesep,setname,filesep,dphfle(r,1).name]);
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
		Db(r) = size(tp_b,1);
		Da = size(tp_a,1);
		
        if ~strcmp(type,'model')
            % calculate model fidelity
            wba = k_b(:,end)./sum(k_b,2);
            fdl(r) = sum(wba)/Db(r)-sum(wba0)/Db0;
            if strcmp(type,'wba')
                simres = ...
                    load([datdir,filesep,setname,filesep,simfle(r,1).name]);
                dt_gt = simres.res.dt_gt;
                for n = 1:numel(dt_gt)
                    dt_gt{n}(:,1) = dt_gt{n}(:,1)*simprm.rate;
                end
                [a0,T0] = calcDPHprm(dt_gt,Da0,Db0);
                plotdat = cat(1,plotdat,bindthist(importdata([datdir,...
                    filesep,setname,filesep,plotfle(r,1).name]),dtbin,...
                    T0{2},a0{2}));
            end
            
        elseif Db(r)==Db0 && Da==Da0 % only for success
            % get computation time
            t_comp(r) = dphres{3}.t_dphtest;
            
            % collect ML-DPH sub-transition matrices (STPMs)
            tau_a = (1./sum(k_a,2))';
            tau_b = (1./sum(k_b,2))';
            STPMs{1} = cat(3,STPMs{1},reorderMat(tp_a,...
                sortStates(repmat(val_a,1,Da),tau_a)));
            STPMs{2} = cat(3,STPMs{2},reorderMat(tp_b,...
                sortStates(repmat(val_b,1,Db(r)),tau_b)));
            
            % collect BW transition probability matrices (TPMs)
            bwfle = [datdir,filesep,filesep,setname,filesep,setname,'_',...
                num2str(r),ext_bw];
            if ~exist(bwfle,'file')
                disp(['No BW results found for set ',setname])
                continue
            end
            bwres = load(bwfle);
            states = dphres{4};
            
            % use same state order (ascending value, then ascending lifetime)
            id_w = sortStates(states,(1./(1-sum(bwres.bwres_w{1},2)))');
            id_wo = sortStates(states,(1./(1-sum(bwres.bwres_wo{1},2)))');
            TPMs = cat(3,TPMs,[simprm.val(id0)',tp0,states(id_w)',...
                reorderMat(bwres.bwres_w{1},id_w),...
                reorderMat(bwres.bwres_w{2}(:,:,1),id_w),...
                reorderMat(bwres.bwres_w{2}(:,:,2),id_w),...
                reorderMat(bwres.bwres_wo{1},id_wo),...
                reorderMat(bwres.bwres_wo{2}(:,:,1),id_wo), ...
                reorderMat(bwres.bwres_wo{2}(:,:,2),id_wo)]);
            BICs = cat(3,BICs,dphres{2});
        end
	end
	
	% calculate std
    switch type
        case 'tau'
            dat = cat(1,dat,[Db0,sort(taub0),mean(Db),std(Db),mean(fdl),...
                std(fdl)]);
        case 'wba'
            dat{1} = cat(1,dat{1},[Db0,mean(wba0),mean(Db),std(Db),...
                mean(fdl),std(fdl)]);
            dat{2} = cat(1,dat{2},...
                {dthiststats(sortrows(plotdat,1),dtbin)});
        case 'model'
            TPMs_std = std(TPMs,[],3);
            TPMs = mean(TPMs,3);
            STPMs = {mean(STPMs{1},3),std(STPMs{1},[],3),...
                mean(STPMs{2},3),std(STPMs{2},[],3)};
            BICs_std = std(BICs,[],3);
            BICs = [mean(BICs,3),BICs_std(:,3:end)];
            dat{1} = {TPMs(:,1),TPMs(:,2:J+1),TPMs(:,J+2),...
                TPMs(:,(J+3):(2*J+2)),TPMs(:,(2*J+3):(3*J+2)),...
                TPMs(:,(3*J+3):(4*J+2)),TPMs(:,(4*J+3):(5*J+2)),...
                TPMs(:,(5*J+3):(6*J+2)),TPMs(:,(6*J+3):(7*J+2)),...
                TPMs_std(:,(J+3):(2*J+2)),TPMs_std(:,(2*J+3):(3*J+2)),...
                TPMs_std(:,(3*J+3):(4*J+2)),TPMs_std(:,(4*J+3):(5*J+2)),...
                TPMs_std(:,(5*J+3):(6*J+2)),TPMs_std(:,(6*J+3):(7*J+2))};
            dat{2} = BICs;
            dat{3} = STPMs;
            dat{4} = {mean(t_comp),std(t_comp)};
            return
            
        case 'none'
            dat = cat(1,dat,[Db0,mean(Db),std(Db),mean(fdl),std(fdl)]);
    end
end
switch type
    case 'none'
        [~,id] = sort(dat(:,1));
        dat = dat(id,:);
    case 'tau'
        [~,id] = sort(dat(:,3));
        dat = dat(id,:);
    case 'wba'
        [~,id] = sort(dat{1}(:,2));
        dat{1} = dat{1}(id,:);
end


function st = dthiststats(dthist,dtbin)
st = [];
maxdt = ceil(max(dthist(:,1))/dtbin);
for bin = 1:maxdt
    id = find(dthist(:,1)==bin*dtbin);
    if ~isempty(id)
        cnt = dthist(id,2)';
        if sum(cnt)>0
            prob = dthist(id,3)';
            prob0 = dthist(id,4)';
            st = cat(1,st,[bin*dtbin,mean(cnt),std(cnt),mean(prob),...
                std(prob),mean(prob0),std(prob0)]);
%              st = cat(1,st,[bin*dtbin,cnt(1),0,prob(1),0,prob0(1),0]);
        end
    end
end


function [a0,T0] = calcDPHprm(dt,Da,Db)
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
prob = zeros(size(dt));
D = size(T,1);
v_e = ones(D,1);
t = v_e-T*v_e;
for n = 1:numel(dt)
    prob(n) = ip*(T^(dt(n)-1))*t;
end


function dthist = bindthist(dthist0,dtbin,T,ip)
dthist = [];
maxdt = ceil(max(dthist0(:,1))/dtbin);
for bin = 1:maxdt
    if bin==1
        lowlim = 0;
    else
        lowlim = bin*dtbin-1;
    end
    id = find(dthist0(:,1)>lowlim & dthist0(:,1)<=bin*dtbin);
    if ~isempty(id)
        prob0 = calcDPHprob(T,ip,dthist0(id,1));
        dthist = cat(1,dthist,[bin*dtbin,sum(dthist0(id,2)),...
            sum(dthist0(id,3)),sum(prob0)]);
    end
end


function mat = reorderMat(mat0,id)
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
% Write data into a .json file using JSON formatting
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
        R = size(dat,1);
        for r = 1:R
            writeJson(f,dat{r,1},ntab+1);
            if r~=R
                fprintf(f,',');
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

