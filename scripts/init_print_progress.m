function [n_exist,n_all,fle,max_l,n_digit] = init_print_progress(src0,...
    f_dat,f_run,d_id)
% [nexist,nall,fle,maxl,ndigit] = init_print_progress(src0,fdat,frun,did)
%
% Initialize log printing.

D = length(f_dat);
fle = repmat("",D,1);
n_exist = zeros(D,1);
n_all = zeros(D,1);
for d = 1:D
    isdat = d_id==d;
    n_all(d) = nnz(isdat);
    n_exist(d) = nnz(cellfun(@(x) exist(x,'file'),f_run(isdat)));
    [src_d,name_d] = fileparts(f_dat{d});
    folder_names = split(src_d((length(src0)+1):end),filesep);
    if length(folder_names)>=2
        fle(d) = string([folder_names{2},filesep,name_d]);
    else
        fle(d) = string(folder_names{1});
    end
end

n_digit = nbdigit(sum(n_all));
max_l = max([strlength(fle); length('TOTAL PROGRESS')]);