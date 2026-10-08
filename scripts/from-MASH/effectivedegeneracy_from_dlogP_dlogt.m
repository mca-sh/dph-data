function [Deff,mindP0,mindP1,mindt1,code] = effectivedegeneracy_from_dlogP_dlogt(varargin)
% [Deff,mindP0,mindP1,mindt1,code] = effectivedegeneracy_from_dlogP_dlogt(Plog,dPlog)
% [Deff,mindP0,mindP1,mindt1,code] = effectivedegeneracy_from_dlogP_dlogt(Plog,dPlog,minprob)
% [Deff,mindP0,mindP1,mindt1,code] = effectivedegeneracy_from_dlogP_dlogt(T,ip,dt)
% [Deff,mindP0,mindP1,mindt1,code] = effectivedegeneracy_from_dlogP_dlogt(T,ip,dt,minprob)
%
% Determines effective state degeneracy encoded in dwell time distribution
% obtained from input DPH model parameters.
% The estimation is based on the number of peaks N_peaks and of inflexion 
% points N_infl found in the log-PMF, using the formula: 
% Deff = 2*N_peaks + N_infl; if N_peaks>0,
% Deff = N_infl + 1;         if N_peaks=0.
% "Fake" inflection points located right bewteen two maxima are ignored.
%
% Plog: [1-by-ndt] or {1-by-V}[1-by-ndtv] already calculated log-binned 
%       dwell time PDF
% T: [D-by-D] or {1-by-V}[Dv-by-Dv] true transition probabilities between 
%    degenrate states 
% t: [D-by-1] or {1-by-V}[Dv-by-1] true exit probabilities of degenerate 
%    states
% ip: [1-by-D] or {1-by-V}[1-by-Dv] true starting probabilities of 
%     degenerate states
% minprob: minimum PMF value to calculate
% Deff: scalar or [1-by-V] effective state degeneracies
% mindP0: [1-by-V] min. gap between two adjacent local extrema in PMF
% mindP1: [1-by-V] min. gap between two adjacent local extrema in first 
%         derivative of PMF
% mindt1: [1-by-V] min. gap between two adjacent local extrema in first 
%         derivative of PMF
% code: {1-by-V}[1-by-C] code for distribution shape (0 fo global max in 
%       PMF, 1 for local max in PMF and 2 fo inflexion in PMF)

% set MATLAB search path
msrc = fileparts(mfilename('fullpath'));
addpath(genpath(msrc));


% CONSTANTS
LOGPMIN = 100*eps;

% collect input
if nargin==2
    Plog = varargin{1};
    dPlog = varargin{2};
    pmin = 0;
    
elseif (nargin==3 && sqrt(numel(varargin{1}))==numel(varargin{2})) || ...
        nargin==4
    T = varargin{1};
    ip = varargin{2};
    dt = varargin{3};
    if nargin==4
        pmin = varargin{4};
    else
        pmin = 0;
    end

    % calculate PDF and log-binned histogram
    if ~iscell(T)
        T = {T};
    end
    if ~iscell(ip)
        ip = {ip};
    end
    V = numel(T);
    Plog = cell(1,V);
    dPlog = cell(1,V);
    for v = 1:V
        [Plog{v},~,dPlog{v}] = calc_DPH_prob(T{v},ip{v},dt);
    end

elseif nargin==3
    Plog = varargin{1};
    dPlog = varargin{2};
    pmin = varargin{3};
    
else
    disp('effectivedegeneracy: 2, 3 or 4 input arguments are required.');
    Deff = 0;
    return
end
if ~iscell(Plog)
    Plog = {Plog};
end
if ~iscell(dPlog)
    dPlog = {dPlog};
end

V = numel(Plog);
Deff = zeros(1,V);
mindP0 = NaN(1,V);
mindP1 = NaN(1,V);
mindt1 = NaN(1,V);
code = cell(1,V);
for v = 1:V
    
    % look for a local extremum in log(P)
    if ~isempty(Plog{v})
        ismin0 = false(size(Plog{v}));
        ismax0 = false(size(Plog{v}));
        ok0 = abs(Plog{v})>LOGPMIN & (10.^Plog{v})>pmin;
        ok0([1:2,(end-1):end]) = false; % exclude edges corrupted by diff
        [ismin0(ok0),ismax0(ok0)] = islocalextr(Plog{v}(ok0));
    else
        ismin0 = false(size(Plog{v}));
        ismax0 = false(size(Plog{v}));
        ok0 = [];
    end

    if ~isempty(dPlog{v})
        % look for a local extremum in d2(log(P))/dt2: inflections
        ok1 = abs(dPlog{v})>LOGPMIN & (10.^Plog{v})>pmin;
        ok1([1:2,(end-1):end]) = false; % exclude edges corrupted by diff
        ismin1 = false(size(dPlog{v}));
        ismax1 = false(size(dPlog{v}));
        [ismin1(ok1),ismax1(ok1)] = islocalextr(dPlog{v}(ok1));
    else
        ismin1 = false(size(Plog{v}));
        ismax1 = false(size(Plog{v}));
    end

    % calculate effective degeneracy using 2 states to describe one
    % maximum and 1 state per inflection
    Deff(v) = 2*sum(ismax0)+sum(ismin1);
    nmax0 = sum(ismax0);
    nmin1 = sum(ismin1);

    if nmax0==0
        % add an initiating state for inflection-only (or no-feature) shapes
        Deff(v) = Deff(v)+1;
        idmax0 = [];
        if nmin1==0
            % break here for no-feature shapes (e. g. single exponential decay)
            continue
        end
    else
        % calculate local/global maxima's amplitudes
        mindP0(v) = Inf;
        idmax0 = find(ismax0);
        idmin0 = find(ismin0);
        for i = idmax0
            if isempty(idmin0) || ~any(idmin0<i)
                mindP0(v) = min(mindP0(v),Plog{v}(i)-Plog{v}(1));
            else
                diffi = (i-idmin0);
                minadj = find(diffi==min(diffi(diffi>0)));
                mindP0(v) = min(mindP0(v),Plog{v}(i)-Plog{v}(idmin0(minadj(1))));
            end
        end
    end

    if nmin1>0
        % calculate durations and amplitudes of inflections
        mindP1(v) = Inf;
        mindt1(v) = Inf;
        idmax1 = find(ismax1);
        idmin1 = find(ismin1);
        for i = idmax1
            if isempty(idmin1) || ~any(idmin1<i)
                mindP1(v) = min(mindP1(v),dPlog{v}(i)-dPlog{v}(1));
                mindt1(v) = min(mindt1(v),i-1);
            else
                diffi = (i-idmin1);
                minadj = find(diffi==min(diffi(diffi>0)));
                mindP1(v) = min(mindP1(v),dPlog{v}(i)-dPlog{v}(idmin1(minadj(1))));
                mindt1(v) = min(mindt1(v),i-idmin1(minadj(1)));
            end
        end
    else
        idmin1 = [];
    end

    % determine shape's code (0 for global max., 1 for local max., 2 for
    % inflection)
    id01 = [idmax0,idmin1];
    codemax0 = ones(1,nmax0); % initialize all max as local
    codemax0(Plog{v}(idmax0)==max(Plog{v}(ok0))) = 0; % find global max.
    codemin1 = 2*ones(1,nmin1); % add inflection
    code01 = [codemax0,codemin1];
    [~,ord] = sort(id01); % sort features in appearing order
    code01 = code01(ord);

    % find and remove "fake" inflection points located bewteen two maxima
    fake_inflection = [false, code01(2:end-1)==2 & ...
        (code01(1:end-2)==0 | code01(1:end-2)==1) & ...
        (code01(3:end)==0 | code01(3:end)==1), false];
    Deff(v) = Deff(v)-nnz(fake_inflection);
    code01(fake_inflection) = [];
    code{v} = code01;
end
if V==1
    code = code{1};
end

% 
% function prob = calcDPHprob(T,t,ip,dt)
% ndt = size(dt,2);
% prob = zeros(1,ndt);
% for n = 1:ndt
%     prob(n) = ip*(T^(dt(n)-1))*t;
% end

