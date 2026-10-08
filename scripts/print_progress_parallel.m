function nb = print_progress_parallel(fle,n_exist,n_all,max_l,n_digit,nb)
% print_progress_parallel(src0,f_dat,f_run,d_id,nb)
%
% Print in command window the progress in parallel iAMM analysis using
% format: [data name]: [nb completed]/[total nb]

D = size(fle,1);
str = repmat("",D,1);
for d = 1:D
    str(d) = string([sprintf(['%',num2str(max_l),'s'],fle(d)),...
        sprintf(': %*i/%*i',n_digit,n_exist(d),n_digit,...
        n_all(d))]);
end

if nb>0
    fprintf(repmat('\b',1,nb));
end
nb1 = fprintf('%s\n',str);
nb2 = fprintf([repmat('-',1,strlength(str(1))),'\n']);
nb3 = fprintf('%*s: %*i/%*i\n\n',max_l,'TOTAL PROGRESS',n_digit,...
    sum(n_exist),n_digit,sum(n_all));
nb = nb1+nb2+nb3;