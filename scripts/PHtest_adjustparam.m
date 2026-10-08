function prm = PHtest_adjustparam(prm, datadir)
if contains(datadir, {'dataset1','dataset2','dataset3'}) % GT: Da=1 Db=1-3
    prm.Dmin = [1,1];
    prm.Dmax = [1,4]; 
elseif contains(datadir, {'dataset4'}) % GT: Da=1 Db=4
    prm.Dmin = [1,1];
    prm.Dmax = [1,5];
elseif contains(datadir, {'EBS-IBS', 'D135'})
    prm.excl = true;
end