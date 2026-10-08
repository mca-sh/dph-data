function [a,eigval] = calcexpweight(Q,pik,tol,PHtype)
%% [a,eigval] = calcexpweight(Q,pi,tol,PHtype)
%
% Calculates parameters in the analytical expression of an:
% - either discrete dwell time PMF p_t=a1*l1^(t-1)+a2*l2^(t-1)+...+an*ln^(t-1)
% - or continuous dwell time PDF f(t)=a1*exp(l1*t)+a2*exp(l2*t)+...an*exp(ln*t)
%
% Q: [K-by-K] generator matrix with rate constants
% pi: [1-by-K] initial state probability
% tol: tolerated deviation between two eigenvalues to be equal
% PHtype: 1 for discrete PH, 2 for continuous PH
% a: [1-by-K] weights a of exponential components
% eigval: [1-by-K] eigenvalues of Q
%%

% intializes output weights
a = [];

% calculates eigenvalues
[~,D,W] = eig(Q);
eigval = diag(D);
eigvalr = rnd2tol(eigval,tol);
        
% check for non-distinct eigenvalues
K = size(Q,1);
if tol~=0 && numel(unique(eigvalr))<K
    eigval = [];
    return
end

% calculates eigenvectors
W = W';
V = inv(W);

% calculates exit probabilities
ratecst = Q;
ratecst(~~eye(K)) = 0;
switch PHtype
    case 1 % DPH
        mu = 1-diag(Q)-sum(ratecst,2); 
    case 2 % CPH
        mu = -diag(Q)-sum(ratecst,2); 
    otherwise
        disp('calcexpweight: unknown distribution type.');
        eigval = [];
        return
end

% calculates weights of exponential components
a = zeros(1,K);
for i = 1:K
    a(i) = sum(mu.*W(i,:)')*sum(pik.*V(:,i)');
end
ar = rnd2tol(a,tol);

% check for null weight
if tol~=0 && any(ar==0)
    a = [];
    eigval = [];
    return
end
