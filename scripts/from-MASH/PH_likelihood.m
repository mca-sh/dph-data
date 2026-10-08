function [logL,xfit,yfit,Ty,Mij] = PH_likelihood(PH,a,T,t,data,mat,Ty,Mij)
% [logL,xfit,yfit,Ty,Mij] = PH_likelihood(PH,a,T,t,data,mat,Ty,Mij)
%
% Calculate the log-likelihood of a PH distribution given a dwell time
% histogram of a state agglomerate.
%
% PH: 1 for discrete PH, 2 for continuous PH
% a: [1-by-D] initial state probabilities
% T: [D-by-D] generator matrix for state agglomerate
% t: [D-by-1] exit probabilities
% data: [2-by-nbins] dwell times and histogram counts
% mat: [2*D-by-2*D] initialized with zeros
% Ty: [D-by-D-by-nbins] initialized with zeros and returned as T^(t-1) for
%     discrete PH, or exp(T*t) for continuous PH
% Mij: [D-by-D-by-nbins] initialized with zeros and returned as the K
%      matrix for discrete PH, or the J matrix for continuous P (cf 
%      Baldt2017)

% intialize output
logL = -Inf;
xfit = [];
yfit = [];

% preliminary calculations of matrices
N = numel(data(1,:));
J = numel(a);
for n = 1:N
    switch PH
        case 1 % DPH
            if data(1,n)>1
                mat(1:J,1:J)= T;
                mat(1:J,(J+1):2*J) = t*a;
                mat((J+1):2*J,1:J) = zeros(J);
                mat((J+1):2*J,(J+1):2*J) = T;
                mat = mat^(data(1,n)-1);
                Ty(:,:,n) = mat(1:J,1:J); % T^y in Baldt2017
                Mij(:,:,n) = mat(1:J,(J+1):end); % matrix K in Baldt2017
            else
                Ty(:,:,n) = T^0;
            end
        case 2 % CPH
            mat(1:J,1:J)= T;
            mat(1:J,(J+1):2*J) = t*a;
            mat((J+1):2*J,1:J) = zeros(J);
            mat((J+1):2*J,(J+1):2*J) = T;
            mat = expm(mat*data(1,n));
            Ty(:,:,n) = mat(1:J,1:J); % exp(T*y) in Baldt2017
            Mij(:,:,n) = mat(1:J,(J+1):end); % matrix J in Baldt2017
    end
end

% calculate likelihood of each histogram point
Ln = zeros(1,N);
for n = 1:N
    switch PH
        case 1 % discrete PH
            Ln(n) = a*(Ty(:,:,n)*t);
            
        case 2 % continuous PH
            Ln(n) = a*(Ty(:,:,n)*t); % Markov jump data
            
%             % time-binned data
%             if data(1,n)==0
%                 continue
%             end
%             CDF2 = 1-a*sum(Ty(:,:,n),2);
%             datexist = data(1,:)==(data(1,n)-1);
%             if any(datexist)
%                 CDF1 = 1-a*sum(Ty(:,:,datexist),2);
%             else
%                 CDF1 = 1-a*sum(expm(T*(data(1,n)-1)),2);
%             end
%             Ln(n) = CDF2-CDF1; 
    end
    if isnan(Ln(n)) || isinf(Ln(n)) || imag(Ln(n))~=0 || Ln(n)<0
        return
    end
end

% calculate summed log-likelihood
datvalid = Ln>0;
xfit = data(1,datvalid);
yfit = Ln(datvalid);
logL = sum(data(2,datvalid).*log(Ln(datvalid)));
if logL==0
    logL = -Inf;
end

