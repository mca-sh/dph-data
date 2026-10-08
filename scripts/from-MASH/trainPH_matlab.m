function [a,T,t,logL,m,nb] = trainPH_matlab(PHtype,a0,T0,P,mute)

% default
M = 1E5; % maximum number of EM iterations
dm_disp = 10; % nb. of iterations between two displays
dL_min = 1E-3; % convergeance criteria on likelihood (faster)
d_min = 1E-8; % convergeance criteria on parameters (from SMACKS)
showparam = false;
plotit = false;

% initialize output
a = [];
T = [];
t = [];
nb = 0;

% pre-allocate memory
J = numel(a0);
nDt = numel(P(1,:));
Ty0 = zeros(J,J,nDt);
Mij0 = zeros(J,J,nDt);
denom = zeros(1,nDt);
mat0 = zeros(2*J,2*J);
Nij0 = zeros(J);
B0 = Nij0(1,:);
Z0 = Nij0(1,:);
Ni0 = Nij0(:,1);
E = eye(J);

% calculate initial likelihood
switch PHtype
    case 1 % DPH
        t0 = 1-sum(T0,2);
    case 2 % CPH
        t0 = -sum(T0,2);
end
% [logL,xfit,yfit] = PH_likelihood(PHtype,a0,T0,t0,P);
[logL,xfit,yfit,Ty,Mij] = PH_likelihood(PHtype,a0,T0,t0,P,mat0,Ty0,Mij0);
if isinf(logL)
    return
end

% plot histogram and inital guess distribution
plotit = plotit & J>1;
if plotit
    fig = figure('windowstyle','docked');
    ax = axes('parent',fig,'xscale','log','yscale','log','nextplot','add');
    scatter(ax,P(1,:),P(2,:)/sum(P(2,:)));
    title(ax,sprintf('D=%i',numel(a0)));
    hp = plot(ax,xfit,yfit,'color','red');
end


N = sum(P(2,:));
m = 0;
m_prev = m;
a = a0;
T = T0;
t = t0;
schm = [0,a0,0;zeros(J,1),T0,t0;zeros(1,J+2)];
schm(~~eye(J+2,J+2)) = 0;
schm = schm>0;
while m<M
    a_prev = a;
    T_prev = T;
    t_prev = t;
    logL_prev = logL;
    
    if plotit
        hp.XData = xfit;
        hp.YData = yfit;
        drawnow;
    end

    % E-step
%     [B,Z,Nij,Ni] = ...
%         PH_Estep(PHtype,a,T,t,P,Ty0,Mij0,mat0,denom,B0,Z0,Ni0,Nij0,E);
    [B,Z,Nij,Ni] = ...
        PH_Estep(PHtype,a,T,t,P,Ty,Mij,denom,B0,Z0,Ni0,Nij0,E);

    % M-step
    [a,T,t] = PH_Mstep(PHtype,B,Z,Nij,Ni,N,schm);

    % likelihood
%     [logL,xfit,yfit] = PH_likelihood(PHtype,a,T,t,P);
    [logL,xfit,yfit,Ty,Mij] = PH_likelihood(PHtype,a,T,t,P,mat0,Ty0,Mij0);
    if isinf(logL)
        a = [];
        T = [];
        t = [];
        logL = -Inf;
        if ~mute
            nb = dispProgress('invalid distribution\n',nb);
        end
        break;
    end

    % check for convergence
    dmax = max([max(max(abs(T-T_prev))),max(abs(a-a_prev)),...
        max(abs(t-t_prev))]);
    dL = (logL-logL_prev)/N;
%     dL = logL-logL_prev;
    if dL<0
        T = [];
        a = [];
        t = [];
        logL = -Inf;
        if ~mute
            nb = dispProgress('likelihood decreases.\n',nb);
        end
        break;
    end
    cvg = dL<dL_min | dmax<d_min;

    % show progress
    if ~mute && (cvg || (m-m_prev)>=dm_disp)
        if showparam
            nb = dispProgress(...
                sprintf(['iteration %i: d=%.3E dL=%.3E\n',...
                repmat('%.4f ',1,J),'\n\n',...
                repmat([repmat('%.4f ',1,J),'\n'],1,J)],m,dmax,dL,a,T'),nb);
        else
            nb = dispProgress(sprintf('iteration %i: d=%.3E dL=%.3E\n',...
                m,dmax,dL),nb);
        end
        m_prev = m;
    end
    if cvg
        break
    end

    m = m+1;
end

if m>=M
    a = [];
    T = [];
    t = [];
    logL = -Inf;
    if ~mute
        nb = dispProgress(['maximum number of iterations has been ',...
            'reached\n'],nb);
    end
end

if plotit
    if ~isinf(logL)
        stophere = true;
    end
    close(fig);
end


% function [B,Z,Nij,Ni] = PH_Estep(...
%     PHtype,a,T,t,data,Ty,Mij,mat,denom,B,Z,Ni,Nij,E)
function [B,Z,Nij,Ni] = PH_Estep(...
    PH_type,a,T,t,data,Ty,Mij,denom,B,Z,Ni,Nij,E)

J = numel(a);
N = numel(data(1,:));

switch PH_type
    case 1 % discrete PH
        % preliminary calculations
        for n = 1:N
%             if data(1,n)>1
%                 mat(1:J,1:J)= T;
%                 mat(1:J,(J+1):2*J) = t*a;
%                 mat((J+1):2*J,1:J) = zeros(J);
%                 mat((J+1):2*J,(J+1):2*J) = T;
%                 mat = mat^(data(1,n)-1);
%                 Ty(:,:,n) = mat(1:J,1:J); % T^y in Baldt2017
%                 Mij(:,:,n) = mat(1:J,(J+1):end); % matrix K in Baldt2017
%             else
%                 Ty(:,:,n) = T^0;
%             end
            denom(n) = a*(Ty(:,:,n)*t);
        end

        % expectation calculation
        for j = 1:J
            for n = 1:N
                if denom(n)<=0
                    continue
                end
                B(j) = B(j)+...
                    data(2,n)*a(j)*(E(j,:)*(Ty(:,:,n)*t))/denom(n);
                for j2 = 1:J
                    if data(1,n)>1
                        Nij(j,j2) = Nij(j,j2)+ ...
                            data(2,n)*(T(j,j2)*Mij(j2,j,n))/denom(n);
                    end
                end
                Ni(j) = Ni(j)+...
                    data(2,n)*a*(Ty(:,:,n)*(E(:,j)*t(j)))/denom(n);
            end
        end
    
    case 2 % continuous PH
        % preliminary calculations
        for n = 1:N
%             mat(1:J,1:J)= T;
%             mat(1:J,(J+1):2*J) = t*a;
%             mat((J+1):2*J,1:J) = zeros(J);
%             mat((J+1):2*J,(J+1):2*J) = T;
%             mat = expm(mat*data(1,n));
%             Ty(:,:,n) = mat(1:J,1:J); % exp(T*y) in Baldt2017
%             Mij(:,:,n) = mat(1:J,(J+1):end); % matrix J in Baldt2017
            denom(n) = a*(Ty(:,:,n)*t);
        end

        % expectation calculation
        for j = 1:J
            for n = 1:N
                if denom(n)<=0
                    continue
                end
                B(j) = B(j)+data(2,n)*a(j)*(E(j,:)*(Ty(:,:,n)*t))/denom(n);
                Z(j) = Z(j)+data(2,n)*Mij(j,j,n)/denom(n);
                for j2 = 1:J
                    if j==j2
                        continue
                    end
                    Nij(j,j2) = Nij(j,j2)+...
                        data(2,n)*T(j,j2)*Mij(j2,j,n)/denom(n);
                end
                Ni(j) = Ni(j)+...
                    data(2,n)*a*(Ty(:,:,n)*(E(:,j)*t(j)))/denom(n);
            end
        end
end


function [a,T,t] = PH_Mstep(PH_type,B,Z,Nij,Ni,totcount,schm)

J = numel(B);
B(~schm(1,2:end-1)) = 0;
a = B/totcount;
a = a/sum(a);

if PH_type==2 % continuous PH
    T = Nij./repmat(Z',[1,J]); 
    t = Ni./Z';
    js = 1:J;
    for j = 1:J
        j2s = js(js~=j);
        T(j,j) = -sum(T(j,j2s))-t(j);
    end
    
    % normalize by one
    tp = [T,t];
    tp(~~eye(J,J+1)) = tp(~~eye(J,J+1))+1;
    tp(~(schm(2:end-1,2:end) | ~~eye(J,J+1))) = 0;
    tp = tp./repmat(sum(tp,2),1,J+1);
    T = tp(:,1:J);
    T(~~eye(J)) = T(~~eye(J))-1;
    t = tp(:,J+1);
    
elseif PH_type==1 % discrete PH
    T = zeros(J);
    t = zeros(J,1);
    for j1 = 1:J
        t(j1) = Ni(j1)/(Ni(j1)+sum(Nij(j1,:)));
        for j2 = 1:J
            T(j1,j2) = Nij(j1,j2)/(Ni(j1)+sum(Nij(j1,:)));
        end
    end
end
