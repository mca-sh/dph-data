function analysismethod = PHtest_getmethodfromcalcmode(calcmode)
% Dtermine which analysis method correspond to input integer code.
%
% calcmode:-2: EM-EMM with positive coefficients
%          -1: EM-EMM with negative coefficients
%           0: test canonical connectivities (+10 for searching for 
%              min. complexity)
%           1: test all connectivities (+10 for searching for min. 
%              complexity)
%           2: test "uncoupled" and "irreversible loop" 
%              connectivities (+10 for searching for min. complexity)
%           3: test "uncoupled" connectivity (+10 for searching for 
%              min. complexity)
%           4: test "coupled"  (+10 for searching for min. 
%              complexity)
%           5: test "uncoupled" and "generalized coxian" connectivity 
%              (+10 for searching for min. complexity)
%           6: test "acyclic" connectivity  (+10 for searching for 
%              min. complexity).
%           7: test "acyclic" connectivity with different state 
%              initiations.
%          10: iEMM from Hines et. al. 2015
%          20: iAMM from Hines et. al. 2015
% analysismethod: 'mlph', 'emexp', 'iemm' or 'iamm'

if calcmode<0
    analysismethod = 'emexp';
elseif calcmode>=10 && calcmode<20
    analysismethod = 'iemm';
elseif calcmode>=20
    analysismethod = 'iamm';
else
    analysismethod = 'mlph';
end