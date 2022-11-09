function DPH_analysisRoutine(rootdir,varargin)
% DPH_analysisRoutine(rootdir)
% DPH_analysisRoutine(rootdir,subfolders)
%
% Perform simulation and analysis routine used in the main article.
%
% rootdir: source directory where data folder are present
% subfolders: data sub-folders to be analyzed

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
if ~isempty(varargin)
    subfolders = varargin{1};
    issubdir = true;
else
    subfolders = {};
    issubdir = false;
end
dircnt = dir(rootdir);
for d = 1:size(dircnt,1)
    if any(contains({'.','..'},dircnt(d,1).name)) || ...
            contains(dircnt(d,1).name,'backup') || ...
            (issubdir && ~contains(subfolders,dircnt(d,1).name))
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


function analyzeDirContent(dircnt,h_fig)

% filter out files and non-data folders
if dircnt.isdir && ~any(contains({'.','..','kinetic model','traces_ASCII'},...
        dircnt.name))
    
    dircnt = dir([dircnt.folder,filesep,dircnt.name]);
    for d = 1:size(dircnt,1)
        analyzeDirContent(dircnt(d,1),h_fig);
    end
    
elseif ~dircnt.isdir
    
    % select .mat preset files missing an associated data folder
    [~,fname,fext] = fileparts(dircnt.name);
    if ~(strcmp(fext,'.mat') && ...
            ~exist([dircnt.folder,filesep,fname],'dir'))
        return
    end
    
    % show action
    disp(['simulate and analyze dataset: ',dircnt.name]);
    
    % project creation and data simulation
    DPHtest_simdata(h_fig,[dircnt.folder,filesep,dircnt.name]);
    
    % determine if BW must be applied 
    if any(contains({'presets_quench','presets_dyn'},fname))
        norate = false;
    else
        norate = true;
    end
    [~,dirname,~] = fileparts(dircnt.name);
    
    % TA analysis of simulated state sequences
    DPHtest_analyzedata(h_fig,[dircnt.folder,filesep,dirname],norate,...
        false);
end
