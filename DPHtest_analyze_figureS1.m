function DPHtest_analyze_figureS1(varargin)
% DPHtest_analyze_figureS1()
%
% Simulate and analyze long trajectories ("fact"-times longer than in the manuscript)

% defaults
N = 100;
srcdir = 'sup-dataset1';
Dmax = 5;
rate = 10;
T = 5;
fact = 10;

% start MASH if necessary
chld = get(0,'children');
h_fig = [];
for c = chld'
    if isprop(c,'Name') && contains(c.Name,'MASH-FRET')
        h_fig = c;
        break
    end
end
if isempty(h_fig)
    h_fig = MASH;
end

% mute actions
h = guidata(h_fig);
muteact = h.mute_actions;
ovask = h.param.OpFiles.overwrite_ask;
ovit = h.param.OpFiles.overwrite;
h.mute_actions = true;
h.param.OpFiles.overwrite_ask = false;
h.param.OpFiles.overwrite = true;
guidata(h_fig,h);

% simulate data and collect dwell times
pname = fileparts(mfilename('fullpath'));
if ~strcmp(pname(end),filesep)
    pname = [pname,filesep];
end
pname = [pname,srcdir,filesep];
flist = dir(pname);
for f = 1:size(flist,1)
    
    % generate dwell times
    [~,rname,fext] = fileparts(flist(f,1).name);
    if ~strcmp(fext,'.mat')
        continue
    end
    pname_mdl = [pname,rname,filesep,'kinetic_model',filesep];
    fname_bic = [rname,'_BIC.txt'];
    fname_fitres = [rname,'_mldphfitres.mat'];
    
    dt = DPHtest_simdt_figureS1(h_fig,[flist(f,1).folder,filesep,...
        flist(f,1).name],rate,N,fact);
    states = unique(dt(:,3));
    V = numel(states);
    for v = 1:V
        dt(dt(:,3)==states(v),3) = v;
    end
    
    % ML-DPH analysis
    [~,mdl,~,~,BIC] = script_findBestModel(dt,Dmax,states,1/rate,1,T,false,...
        pname_mdl);
    
    % export BIC
    if ~exist(pname_mdl,'dir')
        mkdir(pname_mdl)
    end
    fid = fopen([pname_mdl,fname_bic], 'Wt');
    fprintf(fid,['D\t',repmat('BIC(state %0.2f)\t',[1,V]),'\n'],...
        states);
    fprintf(fid,['%i\t',repmat('%d\t',[1,V]),'\n'],[1:Dmax;BIC]);
    fclose(fid);
    
    % export DPH fit parameters and computation time
    save([pname_mdl,fname_fitres],'mdl','-mat')
end

% restore actions
h = guidata(h_fig);
h.mute_actions = muteact;
h.param.OpFiles.overwrite_ask = ovask;
h.param.OpFiles.overwrite = ovit;
guidata(h_fig,h);

disp('analysis for figure S1 completed!')


function dt = DPHtest_simdt_figureS1(h_fig,filepresets,rate,N,fact)
% dt = DPHtest_simdt_figureS1(h_fig,filepresets,rate,N,fact)
%
% Simulate single trajectories from input presets file and collect dwell times of N trajectories.
%
% h_fig: handle to MASH main figure
% filepresets: presets file (*.mat)
% rate: frame rate
% N: number of long trajectories to simulate
% fact: multiplication factor of trajectory lengths
% dt: nDt-by-3 dwell times, molecule indexes, state values
%
% ex: dt = DPHtest_simdt_figureS1(sup-dataset1/presets_I_14_long.mat)

% defaults
nChan = 2;
nL = 1;
Itot = [36,0];
gamma = [1,0];
dE = [0,0];
Bt = [0,0];
isbleach = false;
pixdim = 0.53;
bitrate = 14;
viddim = [256,256];
camnoise = 2;
camprm = [0,0,0,0,0,0; ...
    0,0.067,0.95,0,57.8,0; ...
    0,0,0,0,0,0; ...
    0,0,0,0,0,0; ...
    0,0,0,0,0,0];
bgtype = 1;
bgint = [0,0];

% get interface content
h = guidata(h_fig);

% create new simulation project
pushbutton_newProj_Callback([],1,h_fig);

% set experiment settings
def.nChan = nChan;
def.nL = nL;
def.es = cell(nChan,nL);
[pname,fname,fext] = fileparts(filepresets);
def.es{nChan,nL}.div.projttl = fname;
def.es{nChan,nL}.div.molname = 'unknown';
routinetest_setExperimentSettings(h_fig,def,'sim','');

% import presets and set molecule parameters
h_but = getHandlePanelExpandButton(h.uipanel_S_molecules,h_fig);
if strcmp(h_but.String,char(9660))
    pushbutton_panelCollapse_Callback(h_but,[],h_fig);
end

radiobutton_simCoord_Callback(h.radiobutton_randCoord,[],h_fig);

set(h.checkbox_simPrmFile,'value',1);
checkbox_simPrmFile_Callback(h.checkbox_simPrmFile,{pname,[fname,fext]},...
    h_fig);

set(h.checkbox_photon,'value',1);
checkbox_photon_Callback(h.checkbox_photon,[],h_fig);

set(h.edit_totInt,'string',num2str(Itot(1)));
edit_totInt_Callback(h.edit_totInt,[],h_fig);

set(h.edit_dstrbNoise,'string',num2str(Itot(2)));
edit_dstrbNoise_Callback(h.edit_dstrbNoise,[],h_fig);

set(h.edit_gamma,'string',num2str(gamma(1)));
edit_gamma_Callback(h.edit_gamma,[],h_fig);

set(h.edit_gammaW,'string',num2str(gamma(2)));
edit_gammaW_Callback(h.edit_gammaW,[],h_fig);

set(h.edit_simDeD,'string',num2str(dE(1)));
edit_simDeD_Callback(h.edit_simDeD,[],h_fig);

set(h.edit_simDeA,'string',num2str(dE(2)));
edit_simDeA_Callback(h.edit_simDeA,[],h_fig);

set(h.edit_simBtD,'string',num2str(Bt(1)));
edit_simBtD_Callback(h.edit_simBtD,[],h_fig);

set(h.edit_simBtA,'string',num2str(Bt(2)));
edit_simBtA_Callback(h.edit_simBtA,[],h_fig);

h = guidata(h_fig);
proj = h.param.proj{h.param.curr_proj};
L = fact*2*sum(1./sum(proj.sim.curr.gen_dt{3}{2}.kx(:,:,1),2),1);
t_bleach = L;
set_S_photobleaching(isbleach,t_bleach,h_fig);

% set video parameters
h_but = getHandlePanelExpandButton(h.uipanel_S_videoParameters,h_fig);
if strcmp(h_but.String,char(9660))
    pushbutton_panelCollapse_Callback(h_but,[],h_fig);
end

set(h.edit_length,'string',num2str(L));
edit_length_Callback(h.edit_length,[],h_fig);

set(h.edit_simRate,'string',num2str(rate));
edit_simRate_Callback(h.edit_simRate,[],h_fig);

set(h.edit_pixDim,'string',num2str(pixdim));
edit_pixDim_Callback(h.edit_pixDim,[],h_fig);

set(h.edit_simBitPix,'string',num2str(bitrate));
edit_simBitPix_Callback(h.edit_simBitPix,[],h_fig);

set(h.edit_simMov_w,'string',num2str(viddim(1)));
edit_simMov_w_Callback(h.edit_simMov_w,[],h_fig);

set(h.edit_simMov_h,'string',num2str(viddim(2)));
edit_simMov_h_Callback(h.edit_simMov_h,[],h_fig);

set_S_camNoise(camnoise,camprm,h_fig);

% set experimental setup parameters
h_but = getHandlePanelExpandButton(h.uipanel_S_experimentalSetup,h_fig);
if strcmp(h_but.String,char(9660))
    pushbutton_panelCollapse_Callback(h_but,[],h_fig);
end

set_S_psf(false,[],h_fig);
set_S_defocus(false,[],h_fig);
set_S_spatialBG(bgtype,bgint,[],'',[],h_fig);
set_S_dynamicBG(false,[],h_fig);

% generate data and collect dwell times
dt = [];
nbytes = 0;
for n = 1:N
    
    nbytes = dispProgress(...
        ['process molecule ',num2str(n),'/',num2str(N),'\n'],nbytes);
    
    pushbutton_startSim_Callback(h.pushbutton_startSim,[],h_fig);
    
    h = guidata(h_fig);
    p = h.param;
    proj = p.curr_proj;
    
    seq = p.proj{proj}.sim.prm.res_dat{2}(:,2,1);
    dt_n = getDtFromDiscr(seq,1/rate);
    ndt = size(dt_n,1);
    dt = cat(1,dt,[dt_n(:,1),repmat(n,[ndt,1]),dt_n(:,2)]);
end

% close project
pushbutton_closeProj_Callback([],[],h_fig);

dispProgress(['simulation of ',num2str(size(dt,1)),...
    ' dwell times completed!\n'],nbytes);

