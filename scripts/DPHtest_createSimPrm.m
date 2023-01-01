function DPHtest_createSimPrm(dest,varargin)
% DPHtest_createSimPrm(dest)
% DPHtest_createSimPrm(dest,subfolders)
%
% Create pre-set parameter files used to simulate data in the DPH article.
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
ndtmin = 4; % minimum number of observed dwell times
Db_max = 4; % max. state degeneracy of state b
taub1 = 10; % lifetime of state b1 (in data points)
baseb = 20; % base number used to calculate gaps between b-states' lifetimes (taub(d)=taub(1)*baseb^(d-1))
taub0 = taub1*(baseb.^(0:(Db_max-1))); % lifetimes of b-states
fldr0 = {'dataset1','dataset2','dataset3','dataset4','dataset5',...
        'dataset6','dataset7'};

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
fldr = fldr0{1};
L = Inf;
Da = 1; % state degeneracies of state a
Dbs = 1:Db_max; % state degeneracies of state b
taua = 100; % lifetime of state a (in data points)
wba = 1; % fraction of transitions b-->a over b-->b
for Db = Dbs
    taub = taub0(1:Db); % lifetimes of b-states
    W = getWmat(Da,Db,wba);
    if any(contains(subfolders,fldr))
        createPresetsFile(val,N,L,ndtmin,Da,Db,taua,taub,W,rate,...
            [dest,fldr],sprintf('presets_I_%i%i',Da,Db));
    end
end

% dataset 2 (Da=1, Db=1 to 4, wba=0.8)
fldr = 'dataset2';
wba = 0.8; % fraction of transitions b-->a over b-->b
for Db = Dbs
    taub = taub0(1:Db); % lifetimes of b-states
    W = getWmat(Da,Db,wba);
    if any(contains(subfolders,fldr))
        createPresetsFile(val,N,L,ndtmin,Da,Db,taua,taub,W,rate,...
            [dest,fldr],sprintf('presets_II_%i%i',Da,Db));
    end
end

% dataset 3 (Da=1, Db=3, wba=0 to 1)
fldr = 'dataset3';
Db = 3; % max. state degeneracy of state b
taub = taub0(1:Db); % lifetimes of b-states
wbas = 0:0.1:1; % fraction of transitions b-->a over b-->b
for wba = wbas
    W = getWmat(Da,Db,wba);
    if any(contains(subfolders,fldr))
        createPresetsFile(val,N,L,ndtmin,Da,Db,taua,taub,W,rate,...
            [dest,fldr],strrep(sprintf('presets_II_%.2f',wba),'.',''));
    end
end

% dataset 4 (Da=1, Db=3, wba=1, taub2=taub1 to taub3)
fldr = 'dataset4';
taub2s = logspace(log10(taub(1)),log10(taub(3)),12);
taub2s = taub2s(2:end-1);
wba = 1; % fraction of transitions b-->a over b-->b
W = getWmat(Da,Db,wba);
for taub2 = taub2s
    taub(2) = taub2;
    if any(contains(subfolders,fldr))
        createPresetsFile(val,N,L,ndtmin,Da,Db,taua,taub,W,rate,...
            [dest,fldr],...
            strrep(sprintf('presets_I_%.3f',taub(2)/1000),'.',''));
    end
end

% dataset 5 (Da=1, Db=3, wba=0.8, tau_b2=tau_b1 to tau_b3)
fldr = 'dataset5';
wba = 0.8; % fraction of transitions b-->a over b-->b
W = getWmat(Da,Db,wba);
for taub2 = taub2s
    taub(2) = taub2;
    if any(contains(subfolders,fldr))
        createPresetsFile(val,N,L,ndtmin,Da,Db,taua,taub,W,rate,...
            [dest,fldr],...
            strrep(sprintf('presets_II_%.3f',taub(2)/1000),'.',''));
    end
end

% dataset 6 (Da=2, Db=2, taua=[20,200], taub=[50,500], waa=0)
fldr = 'dataset6';
Da = 2; % state degeneracy of state a
Db = 2; % max. state degeneracy of state b
taua = [20,200];
taub = [50,500];
waa = 0;
L = 2*sum([taua,taub]);
ndtmin = Inf;
W = [0,waa,1-waa,0; waa,0,0,1-waa; 1,0,0,0; 0,1,0,0];
if any(contains(subfolders,fldr))
    createPresetsFile(val,N,L,ndtmin,Da,Db,taua,taub,W,rate,[dest,fldr],...
        'presets_quench');
end

% dataset 7 (Da=2, Db=2, taua=[20,200], taub=[50,500], waa=0.1)
fldr = 'dataset7';
waa = 0.1;
W = [0,waa,1-waa,0; waa,0,0,1-waa; 1,0,0,0; 0,1,0,0];
if any(contains(subfolders,fldr))
    createPresetsFile(val,N,L,ndtmin,Da,Db,taua,taub,W,rate,[dest,fldr],...
        'presets_dyn');
end


function W = getWmat(Da,Db,wba)
% W = getWmat(Da,Db,wba)
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
        W(a,Da+b1) = 1/Db;
        W(Da+b1,a) = wba/Da;
    end
    for b2 = 1:Db
        if b1==b2
            continue
        end
        W(Da+b1,Da+b2) = (1-wba)/(Db-1);
    end
end


function createPresetsFile(stateval,N,L,ndt,Da,Db,taua,taub,W,rate,fldr,fname)
% createPresetsFile(val,N,L,ndt,Da,Db,taua,taub,W,rate,fldr,fname)
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
% fldr: destination folder
% fname: destination file name
%
% example: 
% >> createPresetsFile([0.2,0.8],100,Inf,10,1,3,100,[10,200,4000],...
%     getWmat(1,3,0.8),10,'/homes/github/dph-data','presets_II_13')

if ~exist(fldr,'dir')
    mkdir(fldr);
end
if fldr(end)~=filesep
    fldr = [fldr,filesep];
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
save([fldr,fname,'_simprm.mat'],'N','L','ndt','val','tp','ip','rate',...
    '-mat');

% show success
disp(cat(2,'Presets were successfully written in file: ',...
    [fldr,fname,'_simprm.mat']));
