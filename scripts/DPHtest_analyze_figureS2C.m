function DPHtest_analyze_figureS2C(varargin)
% DPHtest_analyze_figureS2C()
%
% Perform analysis routine used in Figure S2C.

% defaults
srcdir = 'dataset3';
destdir = 'sup-dataset3';

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

% mute all events requiering user action
h = guidata(h_fig);
muteact = h.mute_actions;
ovask = h.param.OpFiles.overwrite_ask;
ovit = h.param.OpFiles.overwrite;
h.mute_actions = true;
h.param.OpFiles.overwrite_ask = false;
h.param.OpFiles.overwrite = true;
guidata(h_fig,h);

% loop analysis through root folder's sub-folders
pname0 = fileparts(mfilename('fullpath'));
if ~strcmp(pname0(end),filesep)
    pname0 = [pname0,filesep];
end
pname = [pname0,srcdir,filesep];
pname_out = [pname0,destdir,filesep];
dircnt = dir(pname);
for d = 1:size(dircnt,1)
    if ~(dircnt(d,1).isdir && ~any(contains({'.','..'},dircnt(d,1).name)))
        continue
    end
    analyzeDirContent(dircnt(d,1),pname_out,h_fig);
end

% restore user actions
h = guidata(h_fig);
h.mute_actions = muteact;
h.param.OpFiles.overwrite_ask = ovask;
h.param.OpFiles.overwrite = ovit;
guidata(h_fig,h);

disp('analysis for figure S2C completed!')


function analyzeDirContent(dircnt,pname_out,h_fig)
pname_in = [dircnt.folder,filesep,dircnt.name];
flist = dir(pname_in);
for f = 1:size(flist,1)
    [~,~,fext] = fileparts(flist(f,1).name);
    if ~strcmp(fext,'.mash')
        continue
    end
    
    pushbutton_openProj_Callback({pname_in,flist(f,1).name},[],h_fig);
    
    destdir = [pname_out,dircnt.name,filesep];
    if ~exist(destdir,'dir')
        mkdir(destdir)
    end
    DPHtest_analyzedata(h_fig,destdir,true,true);
end
