function [schmopt,mdl,schmref_D,schmref_tp] = MLDPH(dt,T,sumexp,Dmax,...
    schmref_D,schmref_tp)

% initializes output
schmopt = [];
mdl = {};
BIC = Inf;

% ML-DPH inference: determine state degeneracy
cvg = false;
D = 1;
while ~cvg && D<=Dmax
    fprintf(['for ', num2str(D),' states:\n']);
    mdl = cat(2,mdl,cell(1,1));
    
    % adjust size of reference schemes table
    if numel(schmref_D)<D
        schmref_D = cat(2,schmref_D,cell(1,D-numel(schmref_D)));
    end
    
    % collect all valid transition schemes
    if sumexp
        ntp = 0;
    else
        ntp = 0:D*(D-1);
    end
    schm_D = [];
    for n = ntp
        
        % from ref table
        if numel(schmref_D{D})>=(n+1) && ~isempty(schmref_D{D}{n+1})
            schm_n = schmref_D{D}{n+1};
        
        % find all possible schemes (can be time consuming)
        else
            [schm_n,schmref_tp] = ...
                getDPHtransSchemesWithNTransProb(n,D,D,schmref_tp);
            if numel(schmref_D{D})<(n+1)
                schmref_D{D} = cat(2,schmref_D{D},...
                    cell(1,(n+1)-numel(schmref_D{D})));
            end
            schmref_D{D}{(n+1)} = schm_n;
        end
        if ~iscell(schm_n) && numel(schm_n)==1 && isinf(schm_n)
            fprintf('no transition scheme\n');
            continue
        end
        if isempty(schm_n)
            fprintf(['maximum model complexity has been reached ',...
                '(Dmax=',num2str(D),'): ML-DPH did not converge\n']);
            break
        end

        % append schemes
        S = numel(schm_n);
        for s = 1:S
            schm_D = cat(3,schm_D,schm_n{s});
        end
    end
    
    % test schemes one by one
    S = size(schm_D,3);
    BIC_D = [];
    for s = 1:S
        fprintf(['scheme ',num2str(s),'/',num2str(S),':\n']);
        mdl_s = script_inferPH(dt,T,schm_D(:,:,s),'');

        % calculate BIC
        nfp = sum(sum(mdl_s.schm))-1;
        mdl_s.BIC = nfp*log(mdl_s.N)-2*mdl_s.logL;

        mdl{D} = cat(1,mdl{D},mdl_s);
        BIC_D = cat(2,BIC_D,mdl_s.BIC);
    end

    % model selection
    BIC_D(isinf(BIC_D)) = Inf;
    [BICmin,sopt] = min(BIC_D);
    if BICmin>BIC
        cvg = true;
        fprintf('ML-DPH successfully converge for D=%i!\n',D-1);
    else
        BIC = BICmin;
        schmopt = mdl{D}(sopt).schm;
        D = D+1;
    end
end

