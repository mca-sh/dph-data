function [dlist,specdat,isspec,dumpname] = PHtest_list_data_folders(src0,...
    inputarg)

% collect names of folders to analyze
specdat = PHtest_getdefdatasetfolders; % default data sub folders
isspec = false;
dumpname = [];
for arg = inputarg
    if iscell(arg{1})
        specdat = arg{1}; % input subfolder
        isspec = true;
    elseif ischar(arg{1})
        dumpname = arg{1};
    end
end

% build a list of analysis folders' contents
dircnt = dir(src0);
dlist = [];
for d = 1:size(dircnt,1)
    if ~dircnt(d,1).isdir || ~contains(dircnt(d,1).name,specdat)
        continue
    end
    dlist = cat(1,dlist,dircnt(d,1));
end