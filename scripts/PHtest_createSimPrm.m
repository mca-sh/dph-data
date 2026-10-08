function def = PHtest_createSimPrm(dest,varargin)
%% PHtest_createSimPrm(dest)
% PHtest_createSimPrm(dest,subfolders)
%
% Create pre-set parameter files used to simulate data in the PH article 
% if it doesn't exist:
% dataset 100 (old dataset1): (Db=1, 1000dt, all connectivities)
% dataset 101: (Db=1, 10000dt, all connectivities)
% dataset 102: (Db=1, 100000dt, all connectivities)
% dataset 103: (Db=1, 1000000dt, all connectivities)
% dataset 200 (old dataset2): (Db=2, 1000dt, all connectivities)
% dataset 201 (old dataset17): (Db=2, 10000dt, all connectivities)
% dataset 202 (old dataset18): (Db=2, 100000dt, all connectivities)
% dataset 203: (Db=2, 1000000dt, all connectivities)
% dataset 300 (old dataset3): (Db=3, 1000dt, all connectivities)
% dataset 301 (old dataset15): (Db=3, 10000dt, all connectivities)
% dataset 302 (old dataset16): (Db=3, 100000dt, all connectivities)
% dataset 303: (Db=3, 1000000dt, all connectivities)
% dataset 400 (old dataset4): (Db=4, 1000dt, random connectivities)
% dataset 401: (Db=4, 10000dt, random connectivities)
% dataset 402: (Db=4, 100000dt, random connectivities)
% dataset 403: (Db=4, 1000000dt, random connectivities)
% dataset 210 (old dataset5): (Db=2, 1000dt, wba=0 to 1, taubi+1=taubi to 10*taubi)
% dataset 211 (old dataset19): (Db=2, 10000dt, wba=0 to 1, taubi+1=taubi to 10*taubi)
% dataset 212 (old dataset20): (Db=2, 100000dt, wba=0 to 1, taubi+1=taubi to 10*taubi)
% dataset 310 (old dataset6): (Db=3, 1000dt, wba=0 to 1, taubi+1=taubi to 10*taubi)
% dataset 311 (old dataset8): (Db=3, 10000dt, wba=0 to 1, taubi+1=taubi to 10*taubi)
% dataset 312 (old dataset9): (Db=3, 100000dt, wba=0 to 1, taubi+1=taubi to 10*taubi)
% dataset 410 (old dataset7): (Db=4, 1000dt, wba=0 to 1, taubi+1=taubi to 10*taubi)
% dataset 411: (Db=4, 10000dt, wba=0 to 1, taubi+1=taubi to 10*taubi)
% dataset 412: (Db=4, 100000dt, wba=0 to 1, taubi+1=taubi to 10*taubi)
% dataset 320 (old dataset10): (Db=3, 1000dt, irreversible loop, taubi+1=taubi to 10*taubi)
% dataset 321 (old dataset11): (Db=3, 10000dt, irreversible loop, taubi+1=taubi to 10*taubi)
% dataset 322 (old dataset12): (Db=3, 100000dt, irreversible loop, taubi+1=taubi to 10*taubi)
% dataset 500 (old dataset13): (Da=2, 1000dt, Db=2, quenched disorder waa=0)
% dataset 600 (old dataset14): (Da=2, 1000dt, Db=2, dynamic disorder waa=0.1)
%
% dest: destination directory
% subfolders: data sub-folders to be processed
%
% example: 
% >> PHtest_createSimPrm('/homes/github/dph-data')

% CONSTANTS
N = 100; % number of trajectories
rate = 10; % frame rate (per second)
val = [0.2,0.7]; % FRET values of state a and b
ndt0 = [10,100,1000,10000]; % minimum number of observed dwell times per trajectory
Dbmax = 4; % max. state degeneracy of state b
taua0 = 100; % lifetime of state a (in data points)
wba0 = 0.8; % fraction of transitions b-->a over b-->b
taub1 = 10; % lifetime of state b1 (in data points)
baseb = 20; % base number used to calculate gaps between b-states' lifetimes (taub(d)=taub(1)*baseb^(d-1))
wbas = 0:0.1:1; % tested fractions of transitions b-->a over b-->b
facttaus = 1:10; % tested multiplication factors for lifetime gap
Nschm = 100; % nb. of random transition schemes
NUM_DT = 1000;
DT_MAX = 1E9;
P_MIN = 1E-9;
reffle = 'from-MASH/ref-table-schemes.mat';
def = struct('N',N,'rate',rate,'val',val,'ndt0',ndt0,'Dbmax',Dbmax,'taua0',...
    taua0,'wba0',wba0,'taub1',taub1,'baseb',baseb,'wbas',wbas,'facttaus',...
    facttaus,'Nschm',Nschm,'reffle',reffle);

% generate dwell time set to calculate the effective degeneracy from
dt_set = unique(logspace(0,log10(DT_MAX),NUM_DT));
    
% calculate lifetimes
taub0 = taub1*(baseb.^(0:(Dbmax)));

% initialize reference schemes
schmD = {};
schm_tp = {};
schm_can = {};
ref.schmD = {};
ref.schm_tp = {};
ref.schm_can = {};

% read reference schemes
mpath = fileparts(mfilename('fullpath'));
reffle = [mpath,filesep,reffle];
if exist(reffle,'file')
    ref = load(reffle);
    if isfield(ref,'schmD')
        schmD = ref.schmD;
    end
    if isfield(ref,'schm_tp')
        schm_tp = ref.schm_tp;
    end
    if isfield(ref,'schm_can')
        schm_can = ref.schm_can;
    end
end

% check existence of directory
if ~exist(dest,'dir')
    disp(['Directory ',dest,' not found: process was aborted.']);
    return
end
if dest(end)~=filesep
    dest = [dest,filesep];
end

% check for sub-folders
if ~isempty(varargin) && ~isempty(varargin{1})
    subfolders = varargin{1};
else
    subfolders = PHtest_getdefdatasetfolders;
end

%% dataset 10x (Da=1, Db=1, all transition schemes)
datset = 'dataset10';
fldrs = subfolders(contains(subfolders,datset));
if ~isempty(fldrs)
    Da = 1;
    Db = 1;
    schm = [];
    for n = 0:Db*(Db-1)
        [schm_n,schm_tp] = PHschm_ntrans(n,Db,Db,schm_tp);
        S = numel(schm_n);
        for s = 1:S
            schm = cat(3,schm,double(schm_n{s}));
        end
    end
    for s = 1:size(schm,3)
        schm_s = schm(:,:,s);
        schm_s = [schm_s(1,1:end-1);...
            [schm_s(2:end-1,end),schm_s(2:(end-1),2:(end-1))]];
        W = getWmat(Da,Db,wba0,schm_s);
        for fldr = fldrs
            fle = [dest,fldr{1},filesep,...
                strrep(sprintf('presets_connex_%.1f',s/10),'.',''),...
                '_simprm.mat'];
            if ~exist(fle,'file')
                m = str2double(fldr{1}(end))+1;
                createPresetsFile(val,N,Inf,ndt0(m),Da,Db,taua0,...
                    taub0(1:Db),W,rate,fle);
            end
        end
    end
end

%% dataset 20x (Da=1, Db=2, all transition schemes)
datset = 'dataset20';
fldrs = subfolders(contains(subfolders,datset));
if ~isempty(fldrs)
    Da = 1;
    Db = 2;
    schm = [];
    for n = 0:Db*(Db-1)
        [schm_n,schm_tp] = PHschm_ntrans(n,Db,Db,schm_tp);
        S = numel(schm_n);
        for s = 1:S
            schm = cat(3,schm,double(schm_n{s}));
        end
    end
    for s = 1:size(schm,3)
        schm_s = schm(:,:,s);
        schm_s = [schm_s(1,1:end-1);...
            [schm_s(2:end-1,end),schm_s(2:(end-1),2:(end-1))]];
        W = getWmat(Da,Db,wba0,schm_s);
        
%         prm.val = val([ones(1,Da),repmat(2,1,Db)]);
%         prm.ip = calcsimip(Da,Db);
%         prm.tp = W./repmat([taua0;taub0(1:Db)'],1,Da+Db);
%         [T,~,ip,~] = collectsdegenparam(prm);
%         Deff = effectivedegeneracy_from_dlogP_dlogt(T,ip,dt_set,P_MIN);
%         if Deff(1)<Da || Deff(2)<Db
%             continue
%         end
        
        for fldr = fldrs
            fle = [dest,fldr{1},filesep,...
                strrep(sprintf('presets_connex_%.1f',s/10),'.',''),...
                '_simprm.mat'];
            if ~exist(fle,'file')
                m = str2double(fldr{1}(end))+1;
                createPresetsFile(val,N,Inf,ndt0(m),Da,Db,taua0,...
                    taub0(1:Db),W,rate,fle);
            end
        end
    end
end

%% dataset 30x (Da=1, Db=3, all transition schemes)
datset = 'dataset30';
fldrs = subfolders(contains(subfolders,datset));
if ~isempty(fldrs)
    Da = 1;
    Db = 3;
    schm = [];
    for n = 0:Db*(Db-1)
        [schm_n,schm_tp] = PHschm_ntrans(n,Db,Db,schm_tp);
        S = numel(schm_n);
        for s = 1:S
            schm = cat(3,schm,double(schm_n{s}));
        end
    end
    for s = 1:size(schm,3)
        schm_s = schm(:,:,s);
        schm_s = [schm_s(1,1:end-1);...
            [schm_s(2:end-1,end),schm_s(2:(end-1),2:(end-1))]];
        W = getWmat(Da,Db,wba0,schm_s);
        
%         prm.val = val([ones(1,Da),repmat(2,1,Db)]);
%         prm.ip = calcsimip(Da,Db);
%         prm.tp = W./repmat([taua0;taub0(1:Db)'],1,Da+Db);
%         [T,~,ip,~] = collectsdegenparam(prm);
%         Deff = effectivedegeneracy_from_dlogP_dlogt(T,ip,dt_set,P_MIN);
%         if Deff(1)<Da || Deff(2)<Db
%             continue
%         end
        
        for fldr = fldrs
            fle = [dest,fldr{1},filesep,...
                strrep(sprintf('presets_connex_%.2f',s/100),'.',''),...
                '_simprm.mat'];
            if ~exist(fle,'file')
                m = str2double(fldr{1}(end))+1;
                createPresetsFile(val,N,Inf,ndt0(m),Da,Db,taua0,...
                    taub0(1:Db),W,rate,fle);
            end
        end
    end
end

%% dataset 40x (Da=1, Db=4, random transition schemes)
datset = 'dataset40';
fldrs = subfolders(contains(subfolders,datset));
if ~isempty(fldrs)
    Da = 1;
    Db = 4;
    schm = [];
    for n = 0:Db*(Db-1)
        [schm_n,schm_tp] = PHschm_ntrans(n,Db,Db,schm_tp);
        S = numel(schm_n);
        for s = 1:S
            schm = cat(3,schm,double(schm_n{s}));
        end
    end

    % append reference file
    if ~isequal(schmD,ref.schmD) || ~isequal(schm_tp,ref.schm_tp) || ...
            ~isequal(schm_can,ref.schm_can)
        save(reffle,'schmD','schm_tp','schm_can','-mat');
    end

    n = 0;
    sid = 1:size(schm,3);
    while n<Nschm
        s = randsample(sid,1);
        sid(sid==s) = [];
        schm_s = schm(:,:,s);
        schm_s = [schm_s(1,1:end-1);...
            [schm_s(2:end-1,end),schm_s(2:(end-1),2:(end-1))]];
        W = getWmat(Da,Db,wba0,schm_s);
        
        prm.val = val([ones(1,Da),repmat(2,1,Db)]);
        prm.ip = calcsimip(Da,Db);
        prm.tp = W./repmat([taua0;taub0(1:Db)'],1,Da+Db);
        [T,~,ip,~] = collectsdegenparam(prm);
        Deff = effectivedegeneracy_from_dlogP_dlogt(T,ip,dt_set,P_MIN);
        if Deff(1)<Da || Deff(2)<Db
            continue
        end
        
        n = n+1;
        for fldr = fldrs
            fle = [dest,fldr{1},filesep,...
                strrep(sprintf('presets_connex_%.2f',n/100),'.',''),...
                '_simprm.mat'];
            if ~exist(fle,'file')
                m = str2double(fldr{1}(end))+1;
                createPresetsFile(val,N,Inf,ndt0(m),Da,Db,taua0,...
                    taub0(1:Db),W,rate,fle);
            end
        end
    end
end

%% dataset 21x (Da=1, Db=2, wba=0 to 1, taub_i+1=taub_i to 10*taub_i)
datset = 'dataset21';
fldrs = subfolders(contains(subfolders,datset));
if ~isempty(fldrs)
    Da = 1; % state degeneracies of state a
    Db = 2; % state degeneracies of state b
    for wba = wbas
        W = getWmat(Da,Db,wba,1-eye(Da+Db));
        for facttau = facttaus
            taub = taub1*(facttau.^(0:Db-1));
            for fldr = fldrs
                fle = [dest,fldr{1},filesep,strrep(sprintf(...
                    'presets_II_%i%i_%.1f_%.2f',Da,Db,facttau/10,wba),...
                    '.',''),'_simprm.mat'];
                if ~exist(fle,'file')
                    m = str2double(fldr{1}(end))+1;
                    createPresetsFile(...
                        val,N,Inf,ndt0(m),Da,Db,taua0,taub,W,rate,fle);
                end
            end
        end
    end
end

%% dataset 31x (Da=1, Db=3, wba=0 to 1, taub_i+1=taub_i to 10*taub_i)
datset = 'dataset31';
fldrs = subfolders(contains(subfolders,datset));
if ~isempty(fldrs)
    Da = 1; % state degeneracies of state a
    Db = 3; % state degeneracies of state b
    for wba = wbas
        W = getWmat(Da,Db,wba,1-eye(Da+Db));
        for facttau = facttaus
            taub = taub1*(facttau.^(0:Db-1));
            for fldr = fldrs
                fle = [dest,fldr{1},filesep,strrep(sprintf(...
                    'presets_II_%i%i_%.1f_%.2f',Da,Db,facttau/10,wba),...
                    '.',''),'_simprm.mat'];
                if ~exist(fle,'file')
                    m = str2double(fldr{1}(end))+1;
                    createPresetsFile(...
                        val,N,Inf,ndt0(m),Da,Db,taua0,taub,W,rate,fle);
                end
            end
        end
    end
end

%% dataset 41x (Da=1, Db=4, wba=0 to 1, taub_i+1=taub_i to 10*taub_i)
datset = 'dataset41';
fldrs = subfolders(contains(subfolders,datset));
if ~isempty(fldrs)
    Da = 1; % state degeneracies of state a
    Db = 4; % state degeneracies of state b
    for wba = wbas
        W = getWmat(Da,Db,wba,1-eye(Da+Db));
        for facttau = facttaus
            taub = taub1*(facttau.^(0:Db-1));
            for fldr = fldrs
                fle = [dest,fldr{1},filesep,strrep(sprintf(...
                    'presets_II_%i%i_%.1f_%.2f',Da,Db,facttau/10,wba),'.',''),...
                    '_simprm.mat'];
                if ~exist(fle,'file')
                    m = str2double(fldr{1}(end))+1;
                    createPresetsFile(...
                        val,N,Inf,ndt0(m),Da,Db,taua0,taub,W,rate,fle);
                end
            end
        end
    end
end

%% dataset 32x (Da=1, Db=3, circular, taub_i+1=taub_i to 10*taub_i)
datset = 'dataset32';
fldrs = subfolders(contains(subfolders,datset));
if ~isempty(fldrs)
    Da = 1; % state degeneracies of state a
    Db = 3; % state degeneracies of state b
    W = circshift(eye(Da+Db),-1,1);
    for facttau = facttaus
        taub = taub1*(facttau.^(0:Db-1));
        for fldr = fldrs
            fle = [dest,fldr{1},filesep,strrep(sprintf(...
                'presets_III_%i%i_%.1f',Da,Db,facttau/10),'.',''),...
                '_simprm.mat'];
            if ~exist(fle,'file')
                m = str2double(fldr{1}(end))+1;
                createPresetsFile(val,N,Inf,ndt0(m),Da,Db,taua0,taub,W,...
                    rate,fle);
            end
        end
    end
end

%% dataset 50x (Da=2, Db=2, taua=[20,200], taub=[50,500], waa=0)
datset = 'dataset50';
fldrs = subfolders(contains(subfolders,datset));
if ~isempty(fldrs)
    Da = 2; % state degeneracy of state a
    Db = 2; % max. state degeneracy of state b
    waa = 0;
    W = [0,waa,1-waa,0; waa,0,0,1-waa; 1,0,0,0; 0,1,0,0];
    for fldr = fldrs
        fle = [dest,fldr{1},filesep,'presets_quench_simprm.mat'];
        if ~exist(fle,'file')
            m = str2double(fldr{1}(end))+1;
            createPresetsFile(val,N,Inf,ndt0(m),Da,Db,[20,200],[50,500],W,...
                rate,fle);
        end
    end
end

%% dataset 60x (Da=2, Db=2, taua=[20,200], taub=[50,500], waa=0.1)
datset = 'dataset60';
fldrs = subfolders(contains(subfolders,datset));
if ~isempty(fldrs)
    waa = 0.1;
    W = [0,waa,1-waa,0; waa,0,0,1-waa; 1,0,0,0; 0,1,0,0];
    for fldr = fldrs
        fle = [dest,fldr{1},filesep,'presets_dyn_simprm.mat'];
        if ~exist(fle,'file')
            createPresetsFile(val,N,Inf,ndt0(m),Da,Db,[20,200],[50,500],W,...
                rate,fle);
        end
    end
end


function W = getWmat(Da,Db,wba,schm)
% W = getWmat(Da,Db,wba,schm)
%
% Calculate normalized transition probabilities for datasets 1 to 5 and 
% return full matrix.
%
% Da: degeneracy of state a
% Db: degeneracy of state b
% wba: fractions of transitions b-->a over all transitions b--> 
%
% example:
% >> W = getWmat(1,3,0.8)
% 
% W =
% 
%          0    0.3333    0.3333    0.3333
%     0.8000         0    0.1000    0.1000
%     0.8000    0.1000         0    0.1000
%     0.8000    0.1000    0.1000         0

J = Da+Db;
W = zeros(J,J);
for b1 = 1:Db
    for a = 1:Da
        if sum(schm(a,(Da+1):(Da+Db)))>0
            W(a,Da+b1) = schm(a,(Da+b1))/sum(schm(a,(Da+1):(Da+Db)));
        end
        if sum(schm(Da+b1,1:Da))>0
            if sum(schm((Da+b1),(Da+1):(Da+Db)))
                W(Da+b1,a) = schm((Da+b1),a)*wba/sum(schm(Da+b1,1:Da));
            else
                W(Da+b1,a) = schm((Da+b1),a)/sum(schm(Da+b1,1:Da));
            end
        end
    end
    for b2 = 1:Db
        if b1==b2
            continue
        end
        if sum(schm(Da+b1,(Da+1):(Da+Db)))>0
            W(Da+b1,Da+b2) = schm(Da+b1,Da+b2)*(1-sum(W(Da+b1,1:Da)))/...
                sum(schm(Da+b1,(Da+1):(Da+Db)));
        end
    end
end


function createPresetsFile(stateval,N,L,ndt,Da,Db,taua,taub,W,rate,fle)
% createPresetsFile(val,N,L,ndt,Da,Db,taua,taub,W,rate,fle)
%
% Export presets to a .mat file.
%
% val: [1-by-J] FRET state values
% N: number of trajectories to simulate
% L: maximum trajectory length
% ndt: minimum number of observed dwell times
% Da: degeneracy of state a
% Db: degeneracy of state b
% taua: [1-by-Da] lifetime of a-states (in data points)
% taub: [1-by-Db] lifetime of b-states (in data points)
% W: [J-by-J] normalized transition probability matrix
% rate: frame rate (per second)
% fle: destination file
%
% example: 
% >> createPresetsFile([0.2,0.8],100,Inf,10,1,3,100,[10,200,4000],...
%     getWmat(1,3,0.8),10,'/homes/github/dph-data/presets_II_13_simprm.mat')

[fldr,~,~] = fileparts(fle);
if ~exist(fldr,'dir')
    mkdir(fldr);
end

% define FRET state values
V = numel(stateval);
val = [];
D = [Da,Db];
for v = 1:V
    val = cat(2,val,repmat(stateval(v),1,D(v)));
end

% initial state probabilities
J = Da+Db;
ip = calcsimip(Da,Db); % initial state probabilities

% define transition rate constants
TAU = repmat([taua';taub'],1,J);
tp = W./TAU;

% export to mat file
save(fle,'N','L','ndt','val','tp','ip','rate','-mat');

% show success
disp(cat(2,'Presets were successfully written in file: ',fle));


function ip = calcsimip(Da,Db)
% Define inital state probabilities in simulated models.
%
% Da: number of states in aggregate a
% Db: number of states in aggregate b

ip = [ones(1,Da),zeros(1,Db)]/Da; 

