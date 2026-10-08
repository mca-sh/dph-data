function [schm,schm_tp] = collectsAllSchm(D,nfp_max,varargin)
% [schm,schm_tp] = collectsAllSchm(D,nfp_max)
% [schm,schm_tp] = collectsAllSchm(D,nfp_max,schm_tp)
% [schm,schm_tp] = collectsAllSchm(D,nfp_max,reffile)
%
% Collects all possible non-redundant transition schemes for an absorbing
% Markov chain between D aggregated states.
%
% D: state aggregate size
% nfp_max: maximum number of free parameters
% schm_tp: {1-by-D_big}{1-by-ntp_big}[D+2-by-D+2-by-Stp] already calculated 
%          transition schemes organized by nb. of transition probabilities.
% reffile: reference .mat file where past colllected schemes were stored
%          (save time), with field:
%   reffile.schm_tp: {1-by-D_big}{1-by-ntp_big}[D+2-by-D+2-by-Stp] already 
%                    calculated transition schemes organized by nb. of 
%                    transition probabilities.
% schm: [D+2-by-D+2-by-S] non-redundant transition schemes

% collects input reference file if any
schm_tp = [];
isrefle = false;
if ~isempty(varargin)
    if iscell(varargin{1})
        schm_tp = varargin{1};
    elseif exist(varargin{1},'file')
        isrefle = true;
        reffle = varargin{1};
        ref = load(reffle);
        if isfield(ref,'schm_tp')
            schm_tp = ref.schm_tp;
        end
    end
end

% collects schemes
schm = [];
for nfp = D:min([nfp_max,(D*(D+1)-1)])
    [schm_nfp,schm_tp] = PHschm_nfp(nfp,false,D,D,schm_tp);
    if ~iscell(schm_nfp)
        schm_nfp = {schm_nfp};
    end
    for s = 1:numel(schm_nfp)
        schm = cat(3,schm,schm_nfp{s});
    end
end

% save schemes to reference file if any
if isrefle
    save(reffle,'schm_tp','-mat');
end
