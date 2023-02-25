function [schm,schm_tp] = getDPHtransSchemesWithNTransProb(ntp,Dmin,Dmax,...
    schm_tp)
% [schm,schm_tp] = getDPHtransSchemesWithNTransProb(ntp,Dmin,Dmax,schm_tp)
%
% Determine all valid transition schemes that use a specific number of free
% parameters.
% Free parameters include strictly positive initial (IP), transition (TP) 
% and exit (EP) probabilities and is calculated as: nfp = (IP-1) + TP + EP
%
% ntp: number transition probabilities between states
% Dmin: minimum number of states
% Dmax: maximum number of states
% schm: {1-by-S} DPH transition schemes
% schm_tp: {1-by-Dbig}{1-by-nbig} already calculated transition schemes
%  for transition probabilities only
% 
% example:
% >> [schm,schm_tp] = getTransSchemesWithNfreePrm(4,2,2,[])

% initialize output
schm = [];

% determine bounds of state degeneracies
if isinf(Dmax)
    disp('getDPHtransSchemesWithNTransProb: infinite bounds.');
    return
end
Dmin = max([Dmin,ceil((sqrt(4*(ntp+1))-1)/2)]);
if Dmax<Dmin
    schm = Inf;
    return
end

% collect already calculated schemes
if numel(schm_tp)<Dmax
    schm_tp = [schm_tp cell(1,Dmax-numel(schm_tp))];
end

% calculate transition schemes
for D = Dmin:Dmax
    [schm_D,schm_tp] = fillDPHtransSchemeWithNTransProb(D,ntp,schm_tp);
    if isempty(schm_D)
        continue
    end
    for s = 1:size(schm_D,3)
        if isempty(schm)
            schm = {};
        end
        schm = cat(2,schm,schm_D(:,:,s));
    end
end


function [schm,schm_tp] = fillDPHtransSchemeWithNTransProb(D,n,schm_tp)
% [schm,schm_tp] = fillDPHtransMatWithNTransProb(D,n,schm_tp)
% 
% Return all valid transition schemes for D degenerate states and n free 
% probabilities
%
% D: state degeneracy
% n: number of transition probabilities
% schm: [(D+2)-by-(D+2)-by-S] transition schemes
% schm_tp: {1-by-Dbig}{1-by-nbig} transition schemes for transition probabilties only

% initializes output
schm = [];

% collect already calculated transition schemes
if numel(schm_tp)<D
    schm_tp = [schm_tp,cell(D-numel(schm_tp))];
end
n0 = numel(schm_tp{D});
if n0<(n+1)
    schm_tp{D} = [schm_tp{D},cell(1,n+1-n0)];
end

% build DPH transition schemes
cmb_all = [];
for n_ip = 1:D 
    
    % add all possible initial state schemes
    ip = getFromTransNbRepart(1,n_ip,D);
    for c_ip = 1:size(ip,1)

        for n_ep = 1:D 

            % get all possible transition schemes
            id = n+1;
            if isempty(schm_tp{D}{id})
                schm_tp{D}{id} = getTransSchemesForNtrans(D,n);
            end
            
            for s = 1:size(schm_tp{D}{id},3)
                
                % check for non-connected states
                if ~all(sum([ip(c_ip,:);schm_tp{D}{id}(:,:,s)],1))
                    continue
                end

                % add all possible exit schemes
                ep = getFromTransNbRepart(1,n_ep,D);
                for c_ep = 1:size(ep,1)

                    % check for absorbing states
                    if ~all(sum([schm_tp{D}{id}(:,:,s),ep(c_ep,:)'],2))
                        continue
                    end
                    
                    % check transition path
                    schm_s = cat(1,[0,ip(c_ip,:),0],...
                        [zeros(D,1),schm_tp{D}{id}(:,:,s),ep(c_ep,:)'],...
                        zeros(1,D+2));
                    if ~isValidTransPath(schm_s)
                        continue
                    end

                    % check if valid state connectivity is already saved
                    duplicate = false;
                    cmb = [sum(schm_s,1);sum(schm_s,2)'];
                    cmb_d = cmb(:,2:end-1);
                    [~,ord2] = sort(cmb_d(2,:));
                    cmb_d = cmb_d(:,ord2');
                    [~,ord1] = sort(cmb_d(1,:));
                    cmb_d = cmb_d(:,ord1');
                    cmb = [cmb(:,1),cmb_d,cmb(:,end)];
                    if ~isempty(cmb_all)
                        duplicate = false;
                        for k = 1:size(cmb_all,3)
                            if all(cmb_all(:,:,k)==cmb,[1,2])
                                duplicate = true;
                                break
                            end
                        end
                    end
                    if duplicate
                        continue 
                    end

                    % save scheme
                    cmb_all = cat(3,cmb_all,cmb);
                    schm = cat(3,schm,schm_s);
                end
            end
        end
    end
end


function schm = getTransSchemesForNtrans(D,n)
        
% get all transition repartition possibilities
rprt = getFromTransNbRepart(D-1,n,D);

% select valid and unique state connectivities
R = size(rprt,1);
cn_all = [];
for to = 1:R
    for from = 1:R

        % state connectivity
        cn = rprt([to,from],:);

        % check if state connectivity is coherent
        if ~isValidConnect(cn)
            continue
        end

        % sort states by transition numbers (1st: to, 2nd: from)
        [~,ord_to] = sort(cn(2,:));
        cn = cn(:,ord_to');
        [~,ord_from] = sort(cn(1,:));
        cn = cn(:,ord_from');

        % check if valid state connectivity is already saved
        duplicate = false;
        if ~isempty(cn_all)
            for k = 1:size(cn_all,3)
                if all(cn_all(:,:,k)==cn,[1,2])
                    duplicate = true;
                    break
                end
            end
        end
        if duplicate
            continue 
        end

        % save state conectivity
        cn_all = cat(3,cn_all,cn);
    end
end

% find corresponding transition matrices
C = size(cn_all,3);
ok = false(1,C);
schm = NaN(D,D,C);
for c = 1:C
    % include no-transition matrix
    if all(all(cn_all(:,:,c)==0))
        schm(:,:,c) = false(D,D);
        ok(c) = true;
        continue
    end

    cn = cn_all(:,:,c);
    pth = false(D,D);
    for j1 = 1:D
        if ~isempty(cn) && all(all(cn==0))
            break
        end
        if cn(1,j1)==0
            continue
        end
        pth_prev = pth;
        cn_prev = cn;
        [pth,cn] = getTransPath(cn,pth,j1);
        if isempty(pth)
            pth = pth_prev;
            cn = cn_prev;
        end
    end
    if ~isempty(cn) && all(all(cn==0))
        schm(:,:,c) = pth;
        ok(c) = true;
    end
end
schm = schm(:,:,ok);


function vects = getFromTransNbRepart(nTrd_max,nTr,D)
% vects = getFromTransNbRepart(nTrd_max,nTr,D)
%
% Returns all possible repartitions of a constant number of transitions 
% from D that can make no more than a maximum number of transition each.
%
% nTrd_max: maximum number of transitions from a single state
% nTr: total number of connexions to partition
% D: number of states
% vects: [nCmb-by-D] possible repartitions of transitions

vects = [];
if nTr<0 || nTrd_max<0 || D<=0
    return
end

% control maximum number of transitions to partition
if (nTr/nTrd_max)>D
    fprintf(['Too many transitions (%i) for %i states with max. %i ',...
        'transitions each.\n'],nTr,D,nTrd_max);
    return
end

for nCxj = min([nTrd_max,nTr]):-1:0
    % control sufficient number of states to partition leftover transitions
    if (nTr-nCxj)>((D-1)*nTrd_max)
        continue
    end
    
    % set transition number for current state
    vect = zeros(1,D);
    vect(1) = nCxj;
    
    % repeat with next states
    if (D-1)>0 && (nTr-nCxj)>=0
        cmb = getFromTransNbRepart(nTrd_max,nTr-nCxj,D-1);
        vect = [repmat(nCxj,size(cmb,1),1),cmb];
    end
    
    % append possible transition repartitions
    vects = cat(1,vects,vect);
end


function ok = isValidConnect(cn0)
% ok = isValidConnect(cn)
%
% Determine whether state connectivities are coherent
%
% cn: [2-by-D] state connectivities (1st: transition from, 2nd: transitions to)
% ok: 1 if state connectivities are valid, 0 otherwise

ok = true;
D = size(cn0,2);

% test state connectivities with circular shift
cn = cn0;
for d1 = 1:D
    % select all possible recieving states
    d2s = circshift(1:D,[0,d1-1]);
    d2s(d2s==d1) = [];
    
    % subtract state 1's connexions to states 2's
    i2 = 0;
    while cn(1,d1)>0 && i2<(D-1)
        i2 = i2+1;
        if cn(2,d2s(i2))>=1
            cn(2,d2s(i2)) = cn(2,d2s(i2))-1;
            cn(1,d1) = cn(1,d1)-1;
        end
    end
end
if all(all(cn==0)) % all connections are fullfilled
    return
end

% test state connectivities without shift
cn = cn0;
for d1 = 1:D
    d2s = 1:D;
    d2s(d2s==d1) = [];
    i2 = 0;
    while cn(1,d1)>0 && i2<(D-1)
        i2 = i2+1;
        if cn(2,d2s(i2))>=1
            cn(2,d2s(i2)) = cn(2,d2s(i2))-1;
            cn(1,d1) = cn(1,d1)-1;
        end
    end
end
if all(all(cn==0)) % all connections are fullfilled
    return
end

% not all state conenctions are fullfilled: incoherent model
ok = false;


function [pth,schm] = getTransPath(schm,pth,d1)

D = size(schm,2);
d2s = [(d1+1):D,1:(d1-1)];
for d2 = d2s
    schm_prev = schm;
    pth_prev = pth;
    if pth(d1,d2)==true % cell in transition matrix already occupied
        schm = schm_prev;
        pth = pth_prev;
        continue
    end
    schm(1,d1) = schm(1,d1)-1;
    schm(2,d2) = schm(2,d2)-1;
    pth(d1,d2) = true;
    if all(all(schm==0))
        break
    end
    if sum(sum(schm<0)) || ...% negative number of remaining transitions
            (sum(schm(1,:)>0)==1 && sum(schm(2,:)>0)==1 && ...
            isequal(schm(1,:),schm(2,:)))% only 1->1 transition possible
        schm = schm_prev;
        pth = pth_prev;
        continue
    end
    
    d1_next = d1+1;
    if d1_next>D
        d1_next = 1;
    end
    while schm(1,d1_next)==0
        d1_next= d1_next+1;
        if d1_next>D
            d1_next = 1;
        end
    end
    [pth,schm] = getTransPath(schm,pth,d1_next);

    if ~isempty(schm) && all(all(schm==0))
        break
    end
    if isempty(schm) || sum(sum(schm<0)) || ...% negative number of remaining transitions
            (sum(schm(1,:)>0)==1 && sum(schm(2,:)>0)==1 && ...
            isequal(schm(1,:),schm(2,:)))% only 1->1 transition possible
        schm = schm_prev;
        pth = pth_prev;
    end
end
if isequal(schm,schm_prev) % dead end
    pth = [];
    schm = [];
end


function ok = isValidTransPath(schm0)
% ok = isValidTransPath(schm)
%
% Check if the transition scheme allows to reach and exit all states.
%
% schm: [D+2-by-D+2] DPH transition scheme
% ok: 1 if scheme is valid, 0 otherwise

% initializes output
ok = false;

D = size(schm0,1)-2;
rched = false(1,D); % 1 if state can be reached, 0 otherwise
exted = false(1,D); % 1 if state can exit, 0 otherwise

% get entry and exit states
d0s = find(schm0(1,2:(end-1)));
exted(~~schm0(2:(end-1),end)') = true;

% follow all possible paths
schm_prev = false(size(schm0));
[rched,exted,~] = followTransPath(0,d0s,schm0,schm_prev,rched,exted);

% check for convergence
if all([rched,exted])
    ok = true;
end


function [rched,exted,schm] = followTransPath(d1,d2s,schm0,schm_prev,rched,...
    exted)

schm = schm_prev;
for d2 = d2s
    
    % mark transition as tested
    schm(d1+1,d2+1) = 1;
    
    % mark states as reached
    D = size(schm0,1)-2;
    d_next = find(schm0(d2+1,2:end));
    rched(d2) = true;
    rched(d_next(d_next<(D+1))) = true;

    % mark states as exited
    if any(d_next==(D+1))
        schm(d2+1,D+2) = 1;
        exted(d2) = true;
        d_next(d_next==(D+1)) = []; % remove absorbing state
    end
    if any(exted(d_next))
        exted(d2) = true;
    end

    % check for convergence
    if all([rched,exted])
        return
    end

    % continue
    if all(schm==schm_prev)
        continue
    end
    schm_prev = schm;

    % follow all possible path
    if isempty(d_next)
        continue
    end
    [rched,exted,schm] = followTransPath(d2,d_next,schm0,schm_prev,rched,...
        exted);
    if d1>0 && any(exted(d_next))
        exted(d1) = true;
    end
end


