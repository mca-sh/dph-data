function DPHtest_createSimPrm(dest,varargin)
% DPHtest_createSimPrm(dest)
% DPHtest_createSimPrm(dest,subfolders)
%
% Create pre-set parameter files used to simulate data in the DPH article 
% if it doesn't exist.
%
% dest: destination directory
% subfolders: data sub-folders to be processed
%
% example: 
% >> DPHtest_createSimPrm('/homes/github/dph-data')

% defaults
N = 100; % number of trajectories
rate = 10; % frame rate (per second)
val = [0.2,0.7]; % FRET values of state a and b
ndtmin = 10; % minimum number of observed dwell times
% ndtmin = Inf; % minimum number of observed dwell times
Db_max = 3; % max. state degeneracy of state b
taub1 = 10; % lifetime of state b1 (in data points)
baseb = 20; % base number used to calculate gaps between b-states' lifetimes (taub(d)=taub(1)*baseb^(d-1))
taub0 = taub1*(baseb.^(0:(Db_max))); % lifetimes of b-states
fldr0 = {'dataset1','dataset2','dataset3','dataset4','dataset5',...
        'dataset6','dataset7','dataset8'};

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
    subfolders = fldr0;
end

% dataset 1 (Da=1, Db=1 to 4, wba=1)
fldr = 'dataset1';
Da = 1; % state degeneracies of state a
Dbs = 1:Db_max; % state degeneracies of state b
taua = 100; % lifetime of state a (in data points)
wba = 1; % fraction of transitions b-->a over b-->b
L = Inf;
for Db = Dbs
    taub = taub0(1:Db); % lifetimes of b-states
%     L = 2*sum([taua,taub]);
    W = getWmat(Da,Db,wba,...
        [[zeros(Da,Da),ones(Da,Db)];[ones(Db,Da),zeros(Db,Db)]]);
    fle = [dest,fldr,filesep,sprintf('presets_I_%i%i_simprm.mat',Da,Db)];
    if any(contains(subfolders,fldr)) && ~exist(fle,'file')
        createPresetsFile(val,N,L,ndtmin,Da,Db,taua,taub,W,rate,fle);
    end
end

% dataset 2 (Da=1, Db=1 to 4, wba=0.5)
fldr = 'dataset2';
wba = 0.8; % fraction of transitions b-->a over b-->b
for Db = Dbs
    taub = taub0(1:Db); % lifetimes of b-states
%     L = 2*sum([taua,taub]);
    W = getWmat(Da,Db,wba,1-eye(Da+Db));
    fle = [dest,fldr,filesep,sprintf('presets_II_%i%i_simprm.mat',Da,Db)];
    if any(contains(subfolders,fldr)) && ~exist(fle,'file')
        createPresetsFile(val,N,L,ndtmin,Da,Db,taua,taub,W,rate,fle);
    end
end

% dataset 3 (Da=1, Db=1 to 4, circular)
fldr = 'dataset3';
for Db = Dbs
    taub = taub0(1:Db); % lifetimes of b-states
%     L = 2*sum([taua,taub]);
    W = circshift(eye(Da+Db),-1,1);
    fle = [dest,fldr,filesep,sprintf('presets_III_%i%i_simprm.mat',Da,Db)];
    if any(contains(subfolders,fldr)) && ~exist(fle,'file')
        createPresetsFile(val,N,L,ndtmin,Da,Db,taua,taub,W,rate,fle);
    end
end

% dataset 4 (Da=1, Db=3, wba=0 to 1)
fldr = 'dataset4';
Db = 2; % max. state degeneracy of state b
taub = taub0(1:Db); % lifetimes of b-states
% L = 2*sum([taua,taub]);
wbas = 0:0.1:1; % fraction of transitions b-->a over b-->b
for wba = wbas
    W = getWmat(Da,Db,wba,1-eye(Da+Db));
    fle = [dest,fldr,filesep,...
        strrep(sprintf('presets_II_%.2f',wba),'.',''),'_simprm.mat'];
    if any(contains(subfolders,fldr)) && ~exist(fle,'file')
        createPresetsFile(val,N,L,ndtmin,Da,Db,taua,taub,W,rate,fle);
    end
end

% dataset 5 (Da=1, Db=3, wba=1, taub2=taub1 to taub3)
fldr = 'dataset5';
taub2s = logspace(log10(taub0(1)),log10(taub0(3)),12);
taub2s = taub2s(2:end-2);
wba = 1; % fraction of transitions b-->a over b-->b
W = getWmat(Da,Db,wba,...
    [[zeros(Da,Da),ones(Da,Db)];[ones(Db,Da),zeros(Db,Db)]]);
for taub2 = taub2s
    taub(2) = taub2;
%     L = 2*sum([taua,taub]);
    fle = [dest,fldr,filesep,...
        strrep(sprintf('presets_I_%.3f',taub(2)/1000),'.',''),...
        '_simprm.mat'];
    if any(contains(subfolders,fldr)) && ~exist(fle,'file')
        createPresetsFile(val,N,L,ndtmin,Da,Db,taua,taub,W,rate,fle);
    end
end

% dataset 6 (Da=1, Db=3, wba=0.8, tau_b2=tau_b1 to tau_b3)
fldr = 'dataset6';
wba = 0.8; % fraction of transitions b-->a over b-->b
W = getWmat(Da,Db,wba,1-eye(Da+Db));
for taub2 = taub2s
    taub(2) = taub2;
%     L = 2*sum([taua,taub]);
    fle = [dest,fldr,filesep,...
        strrep(sprintf('presets_II_%.3f',taub(2)/1000),'.',''),'_simprm.mat'];
    if any(contains(subfolders,fldr)) && ~exist(fle,'file')
        createPresetsFile(val,N,L,ndtmin,Da,Db,taua,taub,W,rate,fle);
    end
end

% dataset 9 (Da=1, Db=2, different transition schemes)
fldr = 'dataset9';
Db = 2;
schmref_tp = [];
schm = [];
for n = 0:2
    [schm_n,schmref_tp] = ...
        getDPHtransSchemesWithNTransProb(n,Db,Db,schmref_tp);
    S = numel(schm_n);
    for s = 1:S
        schm = cat(3,schm,schm_n{s});
    end
end
for s = 1:size(schm,3)
    schm_s = schm(:,:,s);
    schm_s = [schm_s(1,1:end-1);...
        [schm_s(2:end-1,end),schm_s(2:(end-1),2:(end-1))]];
    W = getWmat(Da,Db,wba,schm_s);
    fle = [dest,fldr,filesep,...
        strrep(sprintf('presets_connex_%.1f',s/10),'.',''),'_simprm.mat'];
    if any(contains(subfolders,fldr)) && ~exist(fle,'file')
        createPresetsFile(val,N,L,ndtmin,Da,Db,taua,taub,W,rate,fle);
    end
end

% dataset 7 (Da=2, Db=2, taua=[20,200], taub=[50,500], waa=0)
fldr = 'dataset7';
Da = 2; % state degeneracy of state a
Db = 2; % max. state degeneracy of state b
taua = [20,200];
taub = [50,500];
waa = 0;
L = 2*sum([taua,taub]);
ndtmin = Inf;
W = [0,waa,1-waa,0; waa,0,0,1-waa; 1,0,0,0; 0,1,0,0];
fle = [dest,fldr,filesep,'presets_quench_simprm.mat'];
if any(contains(subfolders,fldr)) && ~exist(fle,'file')
    createPresetsFile(val,N,L,ndtmin,Da,Db,taua,taub,W,rate,fle);
end

% dataset 8 (Da=2, Db=2, taua=[20,200], taub=[50,500], waa=0.1)
fldr = 'dataset8';
waa = 0.1;
W = [0,waa,1-waa,0; waa,0,0,1-waa; 1,0,0,0; 0,1,0,0];
fle = [dest,fldr,filesep,'presets_dyn_simprm.mat'];
if any(contains(subfolders,fldr)) && ~exist(fle,'file')
    createPresetsFile(val,N,L,ndtmin,Da,Db,taua,taub,W,rate,fle);
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
ip = ones(1,J)/J; % initial state probabilities

% define transition rate constants
TAU = repmat([taua';taub'],1,J);
tp = W./TAU;

% export to mat file
save(fle,'N','L','ndt','val','tp','ip','rate','-mat');

% show success
disp(cat(2,'Presets were successfully written in file: ',fle));
