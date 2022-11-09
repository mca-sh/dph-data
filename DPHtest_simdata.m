function DPHtest_simdata(h_fig,filepresets,varargin)
% DPHtest_simdata(h_fig,filepresets)
% DPHtest_simdata(h_fig,filepresets,fact)
%
% Simulate data from input presets file.
%
% h_fig: handle to MASH main figure
% filepresets: presets file (*.mat)
% fact: multiplication factor for the trace length
%
% ex: DPHtest_simdata(dataset2/presets_II_13.mat)

% defaults
nRep = 3;
nChan = 2;
nL = 1;
Itot = [36,0];
gamma = [1,0];
dE = [0,0];
Bt = [0,0];
isbleach = false;
rate = 10;
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
expopt = [false,false,false,true,false,true,false];

if ~isempty(varargin)
    fact = varargin{1};
else
    fact = 1;
end

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

% set export options
h_but = getHandlePanelExpandButton(h.uipanel_S_exportOptions,h_fig);
if strcmp(h_but.String,char(9660))
    pushbutton_panelCollapse_Callback(h_but,[],h_fig);
end

set_S_fileExport(expopt,h_fig);

set(h.popupmenu_opUnits,'value',2);
popupmenu_opUnits_Callback(h.popupmenu_opUnits,[],h_fig);

% generate and export data nRep times
if pname(end)~=filesep
    pname = [pname,filesep];
end
pname_out = [pname,fname];
if ~exist(pname_out,'dir')
    mkdir(pname_out);
end
if pname_out(end)==filesep
    pname_out(end) = [];
end
for r = 1:nRep
    
    if nRep>1
        % set project title
        fname_out = [fname,'_',num2str(r)];
        h = guidata(h_fig);
        p = h.param;
        proj = p.curr_proj;
        p.proj{proj}.exp_parameters{1,2} = fname_out;
    else
        fname_out = fname;
    end

    % generate data
    pushbutton_startSim_Callback(h.pushbutton_startSim,[],h_fig);

    % export trajectories and simulation parameters to ASCII files
    pushbutton_exportSim_Callback({pname_out,fname_out},[],h_fig);

    % save project
    pushbutton_saveProj_Callback({pname_out,[fname_out,'.mash']},[],h_fig);
end


