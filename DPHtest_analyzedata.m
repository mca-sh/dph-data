function DPHtest_analyzedata(h_fig,rootdir,norate,sumexp)
% DPHtest_analyzedata(h_fig,rootdir,norate,sumexp)
%
% Perform DPH data analysis on trajectory files.
%
% h_fig: handle to MASH main figure
% rootdir: dump location
% norate: (1) to omit rate constant inferrence, (0) to infer and export rate constants
% sumexp: (1) to fit sum of exponentials, (0) for DPH
%
% ex: DPHtest_analyzedata(~/Documents/MyDataFolder/DPH-test/dataset2/presets_II_13/traces_ASCII)

% defaults
fretaxis = [-0.2,0.01,1.2];
tdpdiag = 1;
sglcnt = 0;
tdprecalc = 0;
tdpconv = 1;
tdpnorm = 1;
excldt = false;
rearrseq = false;
dmax = 5;
dtbin = 1;
Tdph = 5;
Tbw = 5;

% get interface content
h = guidata(h_fig);

% correct path
if rootdir(end)==filesep
    rootdir(end) = [];
end

% set root folder
pushbutton_rootFolder_Callback({rootdir},[],h_fig);

% set module
switchPan(h.togglebutton_TA,[],h_fig);

% select FRET data and set TDP parameters
h_but = getHandlePanelExpandButton(h.uipanel_TA_transitionDensityPlot,...
    h_fig);
if strcmp(h_but.String,char(9660))
    pushbutton_panelCollapse_Callback(h_but,[],h_fig);
end

set_TA_TDP(3,1,[fretaxis,tdpdiag,sglcnt,tdprecalc,tdpconv,tdpnorm],h_fig);

% set state configuration parameters
h_but = getHandlePanelExpandButton(h.uipanel_TA_stateConfiguration,h_fig);
if strcmp(h_but.String,char(9660))
    pushbutton_panelCollapse_Callback(h_but,[],h_fig);
end

set_TA_stateConfig(2,[2,5,0,0],[1,1,1,1],[],h_fig);

% start TDP clustering
pushbutton_TDPupdateClust_Callback(h.pushbutton_TDPupdateClust,[],h_fig);

% set dwell time histogram parameters
h_but = getHandlePanelExpandButton(h.uipanel_TA_dtHistograms,h_fig);
if strcmp(h_but.String,char(9660))
    pushbutton_panelCollapse_Callback(h_but,[],h_fig);
end

set(h.edit_TA_slBin,'string',num2str(0.01));
edit_TA_slBin_Callback(h.edit_TA_slBin,[],h_fig);

set(h.checkbox_TA_slExcl,'value',excldt);
checkbox_TA_slExcl_Callback(h.checkbox_TA_slExcl,[],h_fig);

set(h.checkbox_tdp_rearrSeq,'value',rearrseq);
checkbox_tdp_rearrSeq_Callback(h.checkbox_tdp_rearrSeq,[],h_fig);

% set kinetic model parameters
h_but = getHandlePanelExpandButton(h.uipanel_TA_kineticModel,h_fig);
if strcmp(h_but.String,char(9660))
    pushbutton_panelCollapse_Callback(h_but,[],h_fig);
end

set_TA_mdl(1,dtbin,dmax,Tdph,Tbw,h_fig);

% start ML-DPH inference
pushbutton_TA_fitMLDPH_Callback({sumexp,true},[],h_fig);

% start BW inference
if ~norate
    pushbutton_TA_refreshModel_Callback(h.pushbutton_TA_refreshModel,...
        h_fig);
end

% save project
[~,dirname,~] = fileparts(rootdir);
pushbutton_saveProj_Callback({rootdir,[dirname,'.mash']},[],h_fig);

% export results
pushbutton_TA_export_Callback(h.pushbutton_TA_export,[],h_fig)
set_TA_expOpt([0,4,0,3,0,0,0,0,0,1,~norate,~norate],h_fig);
pushbutton_expTDPopt_next_Callback({rootdir,[dirname,'.mdl']},[],h_fig);

% close project
pushbutton_closeProj_Callback(h.pushbutton_closeProj,[],h_fig);

