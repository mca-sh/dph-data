function dumpdir = PHtest_getdumpflddir(dumpfldname,datdir)

dumpdir = [];
if isempty(dumpfldname)
    return
end

% sort out dump directories and check for an existing match
cnt = dir(datdir);
isdir = [];
for fld = cnt'
    if any(contains({'..','.'},fld.name))
        isdir = cat(2,isdir,false);
    else
        isdir = cat(2,isdir,fld.isdir);
        if fld.isdir && endsWith(fld.name,dumpfldname)
            dumpdir = fld.name;
            return
        end
    end
end
cnt = cnt(~~isdir);

% build dump directory's number
fldnum = [];
for fld = cnt'
    nameparts = split(fld.name,'-');
    if isnan(str2double(nameparts{1}))
        continue
    end
    fldnum = cat(1,fldnum,nameparts{1});
end
if isempty(fldnum)
    return
end
fldnum = sortrows(fldnum);
dumpdirnum = num2str(str2double(fldnum(end,:))+1);
dumpdir = [repmat('0',length(fldnum(end,:))-length(dumpdirnum)),dumpdirnum,...
    '-',dumpfldname];
