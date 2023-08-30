function ok = cellstrcmp(cellstr,str)
% ok = cellstrcmp(cellstr,str)
%
% Looks for the exact string in each cell.
%
% cellstr: string or {1-by-S} cell array containing strings
% str: string to look for
% ok: [1-by-S] booleans true if found, false otherwise

if ~iscell(cellstr)
    cellstr = {cellstr};
end
S = numel(cellstr);
ok = false(1,S);
for s = 1:S
    if strcmp(cellstr{s},str)
        ok(s) = true;
    end
end
