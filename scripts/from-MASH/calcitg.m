function itg = calcitg(a,eigval,PHtype)
%% itg = calcitg(a,eigval,PHtype)
%
% Calculates integrals of geometric compoenents. 
%
% PHtype: 1 for discrete distribution, 2 for continuous
% a: weights of exponential components in spectral decomposition of PDF
% eigval: exponential constants in spectral decomposition of PDF
% itg: integrals
%%

% initializes output
itg = [];

% ensure correct format
if size(eigval,1)~=size(a,1)
    eigval = eigval';
end

switch PHtype
    case 1 % DPH
        itg = a./(1-eigval);
    case 2 % CPH
        itg = -a./eigval;
    otherwise
        disp('dwelltimeanalysis>calcitg: unknown distribution type.')
        return
end