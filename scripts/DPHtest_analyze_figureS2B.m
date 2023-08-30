function DPHtest_analyze_figureS2B(varargin)
% DPHtest_analyze_figureS2B()
%
% Perform simulation and analysis routine used in Figure S2B.

% defaults
srcdir = 'sup-dataset2';

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
pname = fileparts(mfilename('fullpath'));
if ~strcmp(pname(end),filesep)
    pname = [pname,filesep];
end
pname = [pname,srcdir,filesep];
dircnt = dir(pname);
for d = 1:size(dircnt,1)
    if any(contains({'.','..'},dircnt(d,1).name)) || ...
            ~isempty(strfind(dircnt(d,1).name,'backup'))
        continue
    end
    analyzeDirContent(dircnt(d,1),h_fig);
end

% restore user actions
h = guidata(h_fig);
h.mute_actions = muteact;
h.param.OpFiles.overwrite_ask = ovask;
h.param.OpFiles.overwrite = ovit;
guidata(h_fig,h);

disp('analysis for figure S2B completed!')


function analyzeDirContent(dircnt,h_fig)
if dircnt.isdir
    dircnt = dir([dircnt.folder,filesep,dircnt.name]);
    for d = 1:size(dircnt,1)
        analyzeDirContent(dircnt(d,1),h_fig);
    end
else
    [~,fname,fext] = fileparts(dircnt.name);
    if strcmp(fext,'.mat')
        if exist([dircnt.folder,filesep,fname],'dir')
            return
        end
        DPHtest_simdata(h_fig,[dircnt.folder,filesep,dircnt.name],2);

        [~,dirname,~] = fileparts(dircnt.name);
        DPHtest_analyzedata(h_fig,[dircnt.folder,filesep,dirname],true,...
            false);
    end
end
